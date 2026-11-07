#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Probe what a node offers for measuring the energy that a run takes.

``list`` (the default)
    Every energy and power source of this node: whether it can be read, how
    often it updates, how long a read takes and the smallest step it counts in.

``run -- <command>``
    Runs the command (e.g. the proxy) and reports the energy that every source
    counted over it, with cross-checks: the node against the sum of its parts,
    the devices against the accelerator counters of the node, the energy
    counters against their power integrated over time, and the account of the
    job scheduler where it keeps one.

The sources:

* ``cray-pm``: ``/sys/cray/pm_counters`` of HPE Cray EX nodes (node, CPU,
  memory and accelerator energy and power);
* ``powercap``: the RAPL zones of ``/sys/class/powercap`` (CPU packages, DRAM);
* ``hwmon``: the energy and power sensors of ``/sys/class/hwmon`` (GPU drivers,
  CPU energy drivers, the power meters of Grace);
* ``nvml``, ``amdsmi``, ``rocmsmi``, ``levelzero``: the management libraries of
  NVIDIA, AMD and Intel GPUs, loaded at runtime if they are there.

Run it once per kind of node, as the only job on its node and as the user who
runs the simulations (root may read what users cannot), for example:

    srun -N1 -n1 --exclusive python3 energy_probe.py --json lumi-g.json \\
        run -- ./proxyseissol 100000 200 all

It needs nothing but the standard library of Python 3.6 or later.
"""

import argparse
import ctypes
import glob
import json
import os
import platform
import shutil
import socket
import statistics
import subprocess  # nosec B404: runs the command the user gives, and sstat
import sys
import threading
import time

SCHEMA = "seissol-energy-probe/1"

PCI_VENDORS = {"0x10de": "nvidia", "0x1002": "amd", "0x8086": "intel"}

# the PCI classes of GPUs: VGA, 3D and display controllers, and processing accelerators
GPU_CLASSES = ("0x0300", "0x0302", "0x0380", "0x1200")


class Reading:
    """A value of a source: joules for an energy counter, watts for a power sensor.

    ``timestamp`` is the time of the value as the source gives it, in seconds
    from an origin of its own, or None if the source gives none.
    """

    def __init__(self, value, timestamp=None):
        self.value = value
        self.timestamp = timestamp


class Channel:
    """One quantity of one source that can be read on its own."""

    def __init__(self, name, source, quantity, scope, read, info=None, wrap=None):
        # e.g. "cray-pm/energy", "nvml/0000:3b:00.0"
        self.name = name
        self.source = source
        # "energy" for a counter that only grows, in joules; "power" in watts
        self.quantity = quantity
        # what it measures: node, cpu, memory, accelerator, device, package, ...
        self.scope = scope
        self.read = read
        self.info = info or {}
        # the value at which an energy counter starts over at zero, in joules
        self.wrap = wrap
        self.error = None


def read_text(path):
    with open(path) as handle:
        return handle.read().strip()


def sysfs(root, path):
    return os.path.join(root, path.lstrip("/"))


def pci_address(domain, bus, device, function):
    return "{:04x}:{:02x}:{:02x}.{:x}".format(domain, bus, device, function)


# ---------------------------------------------------------------------------
# sysfs sources
# ---------------------------------------------------------------------------


def parse_cray_pm(text):
    """A file of pm_counters: '<value> <unit> <timestamp> us', or a bare value."""
    parts = text.split()
    value = float(parts[0])
    timestamp = None
    if len(parts) >= 4 and parts[3] == "us":
        timestamp = float(parts[2]) * 1e-6
    return value, (parts[1] if len(parts) > 1 else None), timestamp


def cray_pm_channels(root, notes):
    base = sysfs(root, "/sys/cray/pm_counters")
    if not os.path.isdir(base):
        return [], {}
    channels = []
    meta = {}
    for name in sorted(os.listdir(base)):
        path = os.path.join(base, name)
        if name.endswith("energy") or name.endswith("power"):
            quantity = "energy" if name.endswith("energy") else "power"
            prefix = name[: -len(quantity)].rstrip("_")
            if prefix == "":
                scope = "node"
            elif prefix.startswith("accel"):
                scope = "accelerator"
            else:
                scope = prefix

            def read(path=path):
                value, _, timestamp = parse_cray_pm(read_text(path))
                return Reading(value, timestamp)

            channels.append(
                Channel(
                    "cray-pm/" + name, "cray-pm", quantity, scope, read, {"file": path}
                )
            )
        else:
            try:
                meta[name] = read_text(path)
            except OSError as error:
                meta[name] = "unreadable: {}".format(error.strerror)
    if "freshness" in meta and "raw_scan_hz" in meta:
        notes.append(
            "cray-pm: freshness {} at raw_scan_hz {}".format(
                meta["freshness"], meta["raw_scan_hz"]
            )
        )
    return channels, meta


def powercap_channels(root, notes):
    channels = []
    for zone in sorted(glob.glob(sysfs(root, "/sys/class/powercap/*"))):
        energy = os.path.join(zone, "energy_uj")
        if not os.path.exists(energy):
            continue
        try:
            name = read_text(os.path.join(zone, "name"))
        except OSError:
            name = os.path.basename(zone)
        wrap = None
        try:
            wrap = float(read_text(os.path.join(zone, "max_energy_range_uj"))) * 1e-6
        except (OSError, ValueError):
            pass
        if name.startswith("package"):
            scope = "package"
        elif name == "dram":
            scope = "memory"
        else:
            scope = name

        def read(energy=energy):
            return Reading(float(read_text(energy)) * 1e-6)

        channels.append(
            Channel(
                "powercap/{}/{}".format(os.path.basename(zone), name),
                "powercap",
                "energy",
                scope,
                read,
                {"file": energy},
                wrap,
            )
        )
    return channels


def hwmon_channels(root, notes):
    channels = []
    for hwmon in sorted(glob.glob(sysfs(root, "/sys/class/hwmon/hwmon*"))):
        try:
            driver = read_text(os.path.join(hwmon, "name"))
        except OSError:
            driver = os.path.basename(hwmon)
        device = os.path.realpath(os.path.join(hwmon, "device"))
        pci = None
        for part in reversed(device.split(os.sep)):
            if part.count(":") == 2 and "." in part:
                pci = part
                break
        files = sorted(
            glob.glob(os.path.join(hwmon, "energy*_input"))
            + glob.glob(os.path.join(hwmon, "power*_input"))
            + glob.glob(os.path.join(hwmon, "power*_average"))
        )
        for path in files:
            base = os.path.basename(path)
            sensor = base.split("_")[0]
            quantity = "energy" if base.startswith("energy") else "power"
            label = sensor
            try:
                label = read_text(os.path.join(hwmon, sensor + "_label"))
            except OSError:
                pass

            def read(path=path):
                # energy in microjoules, power in microwatts
                return Reading(float(read_text(path)) * 1e-6)

            info = {"file": path, "driver": driver}
            if pci is not None:
                info["pci"] = pci
            channels.append(
                Channel(
                    "hwmon/{}/{}/{}".format(
                        os.path.basename(hwmon), driver, label or base
                    ),
                    "hwmon",
                    quantity,
                    driver,
                    read,
                    info,
                )
            )
    return channels


def pci_gpus(root):
    """The GPUs on the PCI bus, by address, with their vendor."""
    gpus = {}
    for device in sorted(glob.glob(sysfs(root, "/sys/bus/pci/devices/*"))):
        try:
            cls = read_text(os.path.join(device, "class"))
            vendor = read_text(os.path.join(device, "vendor"))
        except OSError:
            continue
        if cls[:6] in GPU_CLASSES and vendor in PCI_VENDORS:
            gpus[os.path.basename(device)] = PCI_VENDORS[vendor]
    return gpus


# ---------------------------------------------------------------------------
# GPU management libraries
# ---------------------------------------------------------------------------


def load_library(candidates):
    errors = []
    for candidate in candidates:
        try:
            return ctypes.CDLL(candidate), candidate
        except OSError as error:
            errors.append(str(error))
    return None, "; ".join(errors[:2])


def rocm_candidates(name):
    paths = [name]
    roots = [os.environ.get("ROCM_PATH", ""), "/opt/rocm"] + sorted(
        glob.glob("/opt/rocm-*"), reverse=True
    )
    for rootdir in roots:
        if rootdir:
            paths.append(os.path.join(rootdir, "lib", name))
    return paths


def nvml_channels(notes):
    lib, where = load_library(["libnvidia-ml.so.1", "libnvidia-ml.so"])
    if lib is None:
        notes.append("nvml: not found ({})".format(where))
        return []
    u32 = ctypes.c_uint
    handle_type = ctypes.c_void_p

    def fn(name, *argtypes):
        function = getattr(lib, name)
        function.restype = ctypes.c_int
        function.argtypes = list(argtypes)
        return function

    init = fn("nvmlInit_v2")
    count_fn = fn("nvmlDeviceGetCount_v2", ctypes.POINTER(u32))
    by_index = fn("nvmlDeviceGetHandleByIndex_v2", u32, ctypes.POINTER(handle_type))
    get_name = fn("nvmlDeviceGetName", handle_type, ctypes.c_char_p, u32)
    energy_fn = fn(
        "nvmlDeviceGetTotalEnergyConsumption",
        handle_type,
        ctypes.POINTER(ctypes.c_ulonglong),
    )
    power_fn = fn("nvmlDeviceGetPowerUsage", handle_type, ctypes.POINTER(u32))
    limit_fn = fn("nvmlDeviceGetEnforcedPowerLimit", handle_type, ctypes.POINTER(u32))
    clock_fn = fn(
        "nvmlDeviceGetClockInfo", handle_type, ctypes.c_int, ctypes.POINTER(u32)
    )

    class PciInfo(ctypes.Structure):
        _fields_ = [
            ("busIdLegacy", ctypes.c_char * 16),
            ("domain", u32),
            ("bus", u32),
            ("device", u32),
            ("pciDeviceId", u32),
            ("pciSubSystemId", u32),
            ("busId", ctypes.c_char * 32),
        ]

    pci_fn = fn("nvmlDeviceGetPciInfo_v3", handle_type, ctypes.POINTER(PciInfo))

    status = init()
    if status != 0:
        notes.append("nvml: nvmlInit_v2 returned {} ({})".format(status, where))
        return []
    count = u32(0)
    if count_fn(ctypes.byref(count)) != 0:
        notes.append("nvml: no device count")
        return []
    channels = []
    for index in range(count.value):
        handle = handle_type()
        if by_index(index, ctypes.byref(handle)) != 0:
            continue
        name = ctypes.create_string_buffer(96)
        get_name(handle, name, 96)
        pci = PciInfo()
        address = "index{}".format(index)
        if pci_fn(handle, ctypes.byref(pci)) == 0:
            address = pci_address(pci.domain, pci.bus, pci.device, 0)
        info = {"device": name.value.decode(errors="replace"), "pci": address}
        limit = u32(0)
        if limit_fn(handle, ctypes.byref(limit)) == 0:
            info["powerLimit"] = limit.value * 1e-3
        clock = u32(0)
        if clock_fn(handle, 1, ctypes.byref(clock)) == 0:
            info["smClockMHz"] = clock.value

        def read_energy(handle=handle):
            value = ctypes.c_ulonglong(0)
            status = energy_fn(handle, ctypes.byref(value))
            if status != 0:
                raise OSError(
                    "nvmlDeviceGetTotalEnergyConsumption returned {}".format(status)
                )
            return Reading(value.value * 1e-3)

        def read_power(handle=handle):
            value = u32(0)
            status = power_fn(handle, ctypes.byref(value))
            if status != 0:
                raise OSError("nvmlDeviceGetPowerUsage returned {}".format(status))
            return Reading(value.value * 1e-3)

        channels.append(
            Channel(
                "nvml/{}/energy".format(address),
                "nvml",
                "energy",
                "device",
                read_energy,
                info,
            )
        )
        channels.append(
            Channel(
                "nvml/{}/power".format(address),
                "nvml",
                "power",
                "device",
                read_power,
                info,
            )
        )
    return channels


def amdsmi_channels(notes):
    lib, where = load_library(rocm_candidates("libamd_smi.so"))
    if lib is None:
        notes.append("amdsmi: not found ({})".format(where))
        return []
    u32 = ctypes.c_uint32
    u64 = ctypes.c_uint64
    handle_type = ctypes.c_void_p

    def fn(name, *argtypes):
        function = getattr(lib, name)
        function.restype = ctypes.c_int
        function.argtypes = list(argtypes)
        return function

    class PowerCapInfo(ctypes.Structure):
        _fields_ = [
            ("power_cap", u64),
            ("default_power_cap", u64),
            ("dpm_cap", u64),
            ("min_power_cap", u64),
            ("max_power_cap", u64),
            ("reserved", u64 * 3),
        ]

    init = fn("amdsmi_init", u64)
    sockets_fn = fn(
        "amdsmi_get_socket_handles", ctypes.POINTER(u32), ctypes.POINTER(handle_type)
    )
    processors_fn = fn(
        "amdsmi_get_processor_handles",
        handle_type,
        ctypes.POINTER(u32),
        ctypes.POINTER(handle_type),
    )
    bdf_fn = fn("amdsmi_get_gpu_device_bdf", handle_type, ctypes.POINTER(u64))
    energy_fn = fn(
        "amdsmi_get_energy_count",
        handle_type,
        ctypes.POINTER(u64),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(u64),
    )
    cap_fn = fn(
        "amdsmi_get_power_cap_info", handle_type, u32, ctypes.POINTER(PowerCapInfo)
    )

    # AMDSMI_INIT_AMD_GPUS
    status = init(1 << 1)
    if status != 0:
        notes.append("amdsmi: amdsmi_init returned {} ({})".format(status, where))
        return []
    count = u32(0)
    sockets_fn(ctypes.byref(count), None)
    sockets = (handle_type * max(count.value, 1))()
    if sockets_fn(ctypes.byref(count), sockets) != 0:
        notes.append("amdsmi: no sockets")
        return []
    channels = []
    for socket_index in range(count.value):
        processor_count = u32(0)
        processors_fn(sockets[socket_index], ctypes.byref(processor_count), None)
        processors = (handle_type * max(processor_count.value, 1))()
        if (
            processors_fn(
                sockets[socket_index], ctypes.byref(processor_count), processors
            )
            != 0
        ):
            continue
        for index in range(processor_count.value):
            handle = handle_type(processors[index])
            bdf = u64(0)
            address = "socket{}-{}".format(socket_index, index)
            if bdf_fn(handle, ctypes.byref(bdf)) == 0:
                value = bdf.value
                address = pci_address(
                    value >> 16, (value >> 8) & 0xFF, (value >> 3) & 0x1F, value & 0x7
                )
            info = {"pci": address}
            cap = PowerCapInfo()
            if cap_fn(handle, 0, ctypes.byref(cap)) == 0:
                info["powerLimit"] = cap.power_cap * 1e-6

            def read(handle=handle):
                accumulator = u64(0)
                resolution = ctypes.c_float(0)
                timestamp = u64(0)
                status = energy_fn(
                    handle,
                    ctypes.byref(accumulator),
                    ctypes.byref(resolution),
                    ctypes.byref(timestamp),
                )
                if status != 0:
                    raise OSError("amdsmi_get_energy_count returned {}".format(status))
                return Reading(
                    accumulator.value * resolution.value * 1e-6, timestamp.value * 1e-9
                )

            channels.append(
                Channel("amdsmi/" + address, "amdsmi", "energy", "device", read, info)
            )
    return channels


def rocmsmi_channels(notes):
    candidates = []
    for name in ("librocm_smi64.so", "librocm_smi64.so.7", "librocm_smi64.so.6"):
        candidates += rocm_candidates(name)
    lib, where = load_library(candidates)
    if lib is None:
        notes.append("rocmsmi: not found ({})".format(where))
        return []
    u32 = ctypes.c_uint32
    u64 = ctypes.c_uint64

    def fn(name, *argtypes):
        function = getattr(lib, name)
        function.restype = ctypes.c_int
        function.argtypes = list(argtypes)
        return function

    init = fn("rsmi_init", u64)
    count_fn = fn("rsmi_num_monitor_devices", ctypes.POINTER(u32))
    pci_fn = fn("rsmi_dev_pci_id_get", u32, ctypes.POINTER(u64))
    energy_fn = fn(
        "rsmi_dev_energy_count_get",
        u32,
        ctypes.POINTER(u64),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(u64),
    )
    cap_fn = fn("rsmi_dev_power_cap_get", u32, u32, ctypes.POINTER(u64))

    status = init(0)
    if status != 0:
        notes.append("rocmsmi: rsmi_init returned {} ({})".format(status, where))
        return []
    count = u32(0)
    if count_fn(ctypes.byref(count)) != 0:
        return []
    channels = []
    for index in range(count.value):
        bdfid = u64(0)
        address = "index{}".format(index)
        if pci_fn(index, ctypes.byref(bdfid)) == 0:
            value = bdfid.value
            # bits 28 to 31 hold the partition, which is not part of the address
            address = pci_address(
                value >> 32, (value >> 8) & 0xFF, (value >> 3) & 0x1F, value & 0x7
            )
        info = {"pci": address}
        cap = u64(0)
        if cap_fn(index, 0, ctypes.byref(cap)) == 0:
            info["powerLimit"] = cap.value * 1e-6

        def read(index=index):
            accumulator = u64(0)
            resolution = ctypes.c_float(0)
            timestamp = u64(0)
            status = energy_fn(
                index,
                ctypes.byref(accumulator),
                ctypes.byref(resolution),
                ctypes.byref(timestamp),
            )
            if status != 0:
                raise OSError("rsmi_dev_energy_count_get returned {}".format(status))
            return Reading(
                accumulator.value * resolution.value * 1e-6, timestamp.value * 1e-9
            )

        channels.append(
            Channel("rocmsmi/" + address, "rocmsmi", "energy", "device", read, info)
        )
    return channels


def levelzero_channels(notes):
    lib, where = load_library(["libze_loader.so.1", "libze_loader.so"])
    if lib is None:
        notes.append("levelzero: not found ({})".format(where))
        return []
    u32 = ctypes.c_uint32
    handle_type = ctypes.c_void_p

    def fn(name, *argtypes):
        function = getattr(lib, name)
        function.restype = ctypes.c_int
        function.argtypes = list(argtypes)
        return function

    class PciAddress(ctypes.Structure):
        _fields_ = [("domain", u32), ("bus", u32), ("device", u32), ("function", u32)]

    class PciSpeed(ctypes.Structure):
        _fields_ = [
            ("gen", ctypes.c_int32),
            ("width", ctypes.c_int32),
            ("maxBandwidth", ctypes.c_int64),
        ]

    class PciProperties(ctypes.Structure):
        _fields_ = [
            ("stype", ctypes.c_int),
            ("pNext", ctypes.c_void_p),
            ("address", PciAddress),
            ("maxSpeed", PciSpeed),
            ("haveBandwidthCounters", ctypes.c_uint8),
            ("havePacketCounters", ctypes.c_uint8),
            ("haveReplayCounters", ctypes.c_uint8),
        ]

    class PowerProperties(ctypes.Structure):
        _fields_ = [
            ("stype", ctypes.c_int),
            ("pNext", ctypes.c_void_p),
            ("onSubdevice", ctypes.c_uint8),
            ("subdeviceId", u32),
            ("canControl", ctypes.c_uint8),
            ("isEnergyThresholdSupported", ctypes.c_uint8),
            ("defaultLimit", ctypes.c_int32),
            ("minLimit", ctypes.c_int32),
            ("maxLimit", ctypes.c_int32),
        ]

    class EnergyCounter(ctypes.Structure):
        _fields_ = [("energy", ctypes.c_uint64), ("timestamp", ctypes.c_uint64)]

    try:
        status = fn("zesInit", u32)(0)
    except AttributeError:
        # a loader older than zesInit: sysman comes with zeInit then
        os.environ["ZES_ENABLE_SYSMAN"] = "1"
        status = fn("zeInit", u32)(0)
    if status != 0:
        notes.append(
            "levelzero: initialization returned {:#x} ({})".format(status, where)
        )
        return []
    drivers_fn = fn("zesDriverGet", ctypes.POINTER(u32), ctypes.POINTER(handle_type))
    devices_fn = fn(
        "zesDeviceGet", handle_type, ctypes.POINTER(u32), ctypes.POINTER(handle_type)
    )
    pci_fn = fn("zesDevicePciGetProperties", handle_type, ctypes.POINTER(PciProperties))
    domains_fn = fn(
        "zesDeviceEnumPowerDomains",
        handle_type,
        ctypes.POINTER(u32),
        ctypes.POINTER(handle_type),
    )
    properties_fn = fn(
        "zesPowerGetProperties", handle_type, ctypes.POINTER(PowerProperties)
    )
    energy_fn = fn(
        "zesPowerGetEnergyCounter", handle_type, ctypes.POINTER(EnergyCounter)
    )
    card_fn = None
    if hasattr(lib, "zesDeviceGetCardPowerDomain"):
        # the domain of the whole card, which loaders before 1.3 do not have
        card_fn = fn(
            "zesDeviceGetCardPowerDomain", handle_type, ctypes.POINTER(handle_type)
        )

    driver_count = u32(0)
    drivers_fn(ctypes.byref(driver_count), None)
    drivers = (handle_type * max(driver_count.value, 1))()
    drivers_fn(ctypes.byref(driver_count), drivers)
    channels = []
    for driver_index in range(driver_count.value):
        device_count = u32(0)
        devices_fn(drivers[driver_index], ctypes.byref(device_count), None)
        devices = (handle_type * max(device_count.value, 1))()
        devices_fn(drivers[driver_index], ctypes.byref(device_count), devices)
        for device_index in range(device_count.value):
            device = handle_type(devices[device_index])
            properties = PciProperties(stype=0x2)
            address = "driver{}-{}".format(driver_index, device_index)
            if pci_fn(device, ctypes.byref(properties)) == 0:
                pci = properties.address
                address = pci_address(pci.domain, pci.bus, pci.device, pci.function)
            domain_count = u32(0)
            domains_fn(device, ctypes.byref(domain_count), None)
            domains = (handle_type * max(domain_count.value, 1))()
            domains_fn(device, ctypes.byref(domain_count), domains)
            parts = []
            for domain_index in range(domain_count.value):
                domain = handle_type(domains[domain_index])
                power = PowerProperties(stype=0xD)
                part = "domain{}".format(domain_index)
                if properties_fn(domain, ctypes.byref(power)) == 0:
                    part = (
                        "tile{}".format(power.subdeviceId)
                        if power.onSubdevice
                        else "package"
                    )
                parts.append((part, domain))
            card = handle_type()
            if card_fn is not None and card_fn(device, ctypes.byref(card)) == 0:
                parts.append(("card", card))
            for part, domain in parts:

                def read(domain=domain):
                    counter = EnergyCounter()
                    status = energy_fn(domain, ctypes.byref(counter))
                    if status != 0:
                        raise OSError(
                            "zesPowerGetEnergyCounter returned {:#x}".format(status)
                        )
                    return Reading(counter.energy * 1e-6, counter.timestamp * 1e-6)

                channels.append(
                    Channel(
                        "levelzero/{}/{}".format(address, part),
                        "levelzero",
                        "energy",
                        {"package": "device", "card": "card"}.get(part, "tile"),
                        read,
                        {"pci": address},
                    )
                )
    return channels


def discover(root):
    notes = []
    channels, cray_meta = cray_pm_channels(root, notes)
    channels += powercap_channels(root, notes)
    channels += hwmon_channels(root, notes)
    for source in (
        nvml_channels,
        amdsmi_channels,
        rocmsmi_channels,
        levelzero_channels,
    ):
        try:
            channels += source(notes)
        except (AttributeError, OSError, ValueError) as error:
            notes.append("{}: {}".format(source.__name__.split("_")[0], error))
    return channels, notes, cray_meta


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------


def try_read(channel):
    try:
        return channel.read()
    except (OSError, ValueError, IndexError) as error:
        channel.error = getattr(error, "strerror", None) or str(error)
        return None


def latency(channel, reads):
    times = []
    for _ in range(reads):
        start = time.perf_counter()
        channel.read()
        times.append(time.perf_counter() - start)
    return statistics.median(times)


def watch(channels, seconds, pause):
    """Polls all channels together; returns per channel the times its value changed."""
    changes = {channel.name: [] for channel in channels}
    steps = {channel.name: [] for channel in channels}
    last = {}
    start = time.perf_counter()
    first = {}
    end_values = {}
    while time.perf_counter() - start < seconds:
        for channel in channels:
            try:
                reading = channel.read()
            except (OSError, ValueError, IndexError):
                continue
            now = time.perf_counter()
            if channel.name not in first:
                first[channel.name] = (now, reading)
            previous = last.get(channel.name)
            if previous is not None and reading.value != previous.value:
                changes[channel.name].append(now)
                if reading.value > previous.value:
                    steps[channel.name].append(reading.value - previous.value)
            last[channel.name] = reading
            end_values[channel.name] = (now, reading)
        time.sleep(pause)
    return changes, steps, first, end_values


def characterize(channels, seconds, pause, reads):
    readable = []
    for channel in channels:
        if try_read(channel) is not None:
            readable.append(channel)
    for channel in readable:
        try:
            channel.info["readLatency"] = latency(channel, reads)
        except (OSError, ValueError, IndexError) as error:
            channel.error = str(error)
    changes, steps, first, end = watch(readable, seconds, pause)
    for channel in readable:
        times = changes[channel.name]
        intervals = [b - a for a, b in zip(times, times[1:])]
        channel.info["changes"] = len(times)
        if intervals:
            channel.info["updateInterval"] = statistics.median(intervals)
        if steps[channel.name]:
            channel.info["smallestStep"] = min(steps[channel.name])
        if (
            channel.quantity == "energy"
            and channel.name in first
            and channel.name in end
        ):
            (t0, r0), (t1, r1) = first[channel.name], end[channel.name]
            elapsed = t1 - t0
            if (
                r0.timestamp is not None
                and r1.timestamp is not None
                and r1.timestamp > r0.timestamp
            ):
                elapsed = r1.timestamp - r0.timestamp
            delta = r1.value - r0.value
            if elapsed > 0 and delta >= 0:
                channel.info["idlePower"] = delta / elapsed
        elif channel.quantity == "power" and channel.name in end:
            channel.info["idlePower"] = end[channel.name][1].value
    return readable


class Sampler(threading.Thread):
    """Reads the given channels every period, for integrating power over time."""

    def __init__(self, channels, period):
        super().__init__(daemon=True)
        self.channels = channels
        self.period = period
        self.samples = {channel.name: [] for channel in channels}
        self.stopped = threading.Event()

    def run(self):
        while not self.stopped.is_set():
            for channel in self.channels:
                try:
                    reading = channel.read()
                except (OSError, ValueError, IndexError):
                    continue
                self.samples[channel.name].append((time.perf_counter(), reading.value))
            self.stopped.wait(self.period)

    def stop(self):
        self.stopped.set()
        self.join()


def integrate(samples):
    """The trapezoidal integral of power samples, in joules, and their time span."""
    total = 0.0
    for (t0, p0), (t1, p1) in zip(samples, samples[1:]):
        total += 0.5 * (p0 + p1) * (t1 - t0)
    span = samples[-1][0] - samples[0][0] if len(samples) > 1 else 0.0
    return total, span


def slurm_step_energy():
    """The energy that Slurm accounts the current job step so far, in joules, or None."""
    job = os.environ.get("SLURM_JOB_ID")
    step = os.environ.get("SLURM_STEP_ID")
    sstat = shutil.which("sstat")
    if not job or step is None or sstat is None:
        return None
    try:
        result = subprocess.run(  # nosec B603: fixed arguments
            [
                sstat,
                "-n",
                "-P",
                "-j",
                "{}.{}".format(job, step),
                "-o",
                "ConsumedEnergyRaw",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    values = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    try:
        return sum(float(value) for value in values) if values else None
    except ValueError:
        return None


def energy_delta(channel, before, after):
    delta = after.value - before.value
    if delta < 0 and channel.wrap:
        delta += channel.wrap
    elapsed = None
    if before.timestamp is not None and after.timestamp is not None:
        elapsed = after.timestamp - before.timestamp
    return delta, elapsed


def run_workload(channels, command, period):
    energy = [channel for channel in channels if channel.quantity == "energy"]
    sampled = list(channels)
    before = {channel.name: try_read(channel) for channel in energy}
    slurm_before = slurm_step_energy()
    sampler = Sampler(sampled, period)
    sampler.start()
    start = time.perf_counter()
    returncode = subprocess.call(
        command
    )  # nosec B603: the command the user asked to run
    wall = time.perf_counter() - start
    sampler.stop()
    after = {channel.name: try_read(channel) for channel in energy}
    slurm_after = slurm_step_energy()

    results = {
        "command": command,
        "returncode": returncode,
        "wallTime": wall,
        "channels": {},
    }
    for channel in channels:
        entry = {}
        if channel.quantity == "energy":
            b, a = before.get(channel.name), after.get(channel.name)
            if b is None or a is None:
                continue
            delta, elapsed = energy_delta(channel, b, a)
            entry["energy"] = delta
            entry["meanPower"] = delta / (elapsed if elapsed else wall)
            if elapsed:
                entry["sensorTime"] = elapsed
        samples = sampler.samples.get(channel.name, [])
        if len(samples) > 1:
            if channel.quantity == "power":
                integral, span = integrate(samples)
                entry["energy"] = integral
                entry["meanPower"] = integral / span if span > 0 else None
            entry["samples"] = len(samples)
        if entry:
            results["channels"][channel.name] = entry
    if slurm_before is not None and slurm_after is not None:
        results["slurmStepEnergy"] = slurm_after - slurm_before
    results["checks"] = cross_checks(
        channels, results["channels"], results.get("slurmStepEnergy")
    )
    return results


def cross_checks(channels, measured, slurm):
    checks = []
    by_name = {channel.name: channel for channel in channels}

    def energy_of(predicate):
        values = [
            entry["energy"]
            for name, entry in measured.items()
            if "energy" in entry
            and by_name[name].quantity == "energy"
            and predicate(by_name[name])
        ]
        return sum(values) if values else None

    node = energy_of(lambda c: c.source == "cray-pm" and c.scope == "node")
    parts = energy_of(lambda c: c.source == "cray-pm" and c.scope not in ("node",))
    accelerators = energy_of(
        lambda c: c.source == "cray-pm" and c.scope == "accelerator"
    )
    if node and parts is not None:
        checks.append(
            "cray-pm: CPU, memory and accelerators are {:.1%} of the node; "
            "the rest is the board, the network card and the fans".format(parts / node)
        )
    for source in ("nvml", "amdsmi", "rocmsmi", "levelzero"):
        devices = energy_of(lambda c, s=source: c.source == s and c.scope == "device")
        if devices is None:
            continue
        if accelerators:
            checks.append(
                "{}: the devices count {:.1%} of the accelerator energy of cray-pm".format(
                    source, devices / accelerators
                )
            )
        zero = [
            name
            for name, entry in measured.items()
            if by_name[name].source == source
            and by_name[name].scope == "device"
            and entry.get("energy", 1.0) <= 0.0
        ]
        if zero:
            checks.append(
                "{}: {} counted no energy, perhaps because another device of the same "
                "package reports for it: {}".format(source, len(zero), ", ".join(zero))
            )
    for name, entry in measured.items():
        channel = by_name[name]
        if channel.source == "nvml" and channel.quantity == "energy":
            power = measured.get(name.replace("/energy", "/power"), {})
            if power.get("energy") and entry.get("energy"):
                checks.append(
                    "{}: energy counter over integrated power {:.3f}".format(
                        name, entry["energy"] / power["energy"]
                    )
                )
    if slurm is not None and node:
        checks.append(
            "slurm: the step accounts {:.1f} J, cray-pm counts {:.1f} J on this node".format(
                slurm, node
            )
        )
    return checks


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------


def host_info(root):
    info = {
        "hostname": socket.gethostname(),
        "kernel": platform.release(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "uid": os.getuid() if hasattr(os, "getuid") else None,
        "cpus": os.cpu_count(),
        "gpus": pci_gpus(root),
        "slurm": {
            key: value
            for key, value in os.environ.items()
            if key
            in (
                "SLURM_JOB_ID",
                "SLURM_STEP_ID",
                "SLURM_CLUSTER_NAME",
                "SLURM_JOB_PARTITION",
            )
        },
    }
    try:
        for line in read_text(sysfs(root, "/proc/cpuinfo")).splitlines():
            if line.startswith("model name") or line.startswith("Model"):
                info["cpu"] = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return info


def format_quantity(value, unit):
    if value is None:
        return "-"
    for factor, prefix in (
        (1e9, "G"),
        (1e6, "M"),
        (1e3, "k"),
        (1, ""),
        (1e-3, "m"),
        (1e-6, "u"),
    ):
        if abs(value) >= factor:
            return "{:.3g} {}{}".format(value / factor, prefix, unit)
    return "{:.3g} {}".format(value, unit)


def print_list(channels, readable, notes):
    names = {channel.name for channel in readable}
    print(
        "{:<52} {:>6} {:>10} {:>10} {:>10} {:>10}".format(
            "source", "kind", "latency", "update", "step", "power"
        )
    )
    for channel in channels:
        if channel.name not in names:
            print(
                "{:<52} {:>6} unreadable: {}".format(
                    channel.name, channel.quantity, channel.error
                )
            )
            continue
        info = channel.info
        unit = "J" if channel.quantity == "energy" else "W"
        print(
            "{:<52} {:>6} {:>10} {:>10} {:>10} {:>10}".format(
                channel.name,
                channel.quantity,
                format_quantity(info.get("readLatency"), "s"),
                format_quantity(info.get("updateInterval"), "s"),
                format_quantity(info.get("smallestStep"), unit),
                format_quantity(info.get("idlePower"), "W"),
            )
        )
    for channel in readable:
        interval = channel.info.get("updateInterval")
        if channel.quantity == "energy" and interval:
            print(
                "hint: {}: epochs of {:.1f} s or more keep its update steps under 5% "
                "of what it counts in an epoch".format(channel.name, 20 * interval)
            )
    for note in notes:
        print("note: " + note)
    if hasattr(os, "getuid") and os.getuid() == 0:
        print(
            "note: running as root, which may read sources that the users of the machine cannot"
        )


def print_run(results):
    print(
        "ran {} in {:.2f} s (exit code {})".format(
            " ".join(results["command"]), results["wallTime"], results["returncode"]
        )
    )
    print("{:<52} {:>12} {:>12}".format("source", "energy", "mean power"))
    for name, entry in sorted(results["channels"].items()):
        print(
            "{:<52} {:>12} {:>12}".format(
                name,
                format_quantity(entry.get("energy"), "J"),
                format_quantity(entry.get("meanPower"), "W"),
            )
        )
    if "slurmStepEnergy" in results:
        print("slurm step: {}".format(format_quantity(results["slurmStepEnergy"], "J")))
    for check in results["checks"]:
        print("check: " + check)


def channel_json(channel, readable):
    entry = {
        "source": channel.source,
        "quantity": channel.quantity,
        "scope": channel.scope,
        "readable": channel in readable,
    }
    if channel.error:
        entry["error"] = channel.error
    if channel.wrap:
        entry["wrap"] = channel.wrap
    entry.update(channel.info)
    return entry


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--json", help="also write the results to this JSON file")
    parser.add_argument(
        "--watch", type=float, default=3.0, help="seconds to poll each source for (3)"
    )
    parser.add_argument(
        "--pause", type=float, default=0.0005, help="seconds between two polls (0.0005)"
    )
    parser.add_argument(
        "--reads", type=int, default=30, help="reads to time per source (30)"
    )
    parser.add_argument(
        "--period",
        type=float,
        default=0.05,
        help="seconds between samples during run (0.05)",
    )
    parser.add_argument(
        "--root", default="/", help="where /sys and /proc are, for testing (/)"
    )
    parser.add_argument("mode", nargs="?", choices=("list", "run"), default="list")
    parser.add_argument(
        "command", nargs=argparse.REMAINDER, help="for run: the command, after --"
    )
    args = parser.parse_args(argv)

    channels, notes, cray_meta = discover(args.root)
    readable = characterize(channels, args.watch, args.pause, args.reads)
    print_list(channels, readable, notes)

    document = {
        "schema": SCHEMA,
        "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": host_info(args.root),
        "crayPm": cray_meta,
        "channels": {
            channel.name: channel_json(channel, readable) for channel in channels
        },
        "notes": notes,
    }
    if args.mode == "run":
        command = args.command[1:] if args.command[:1] == ["--"] else args.command
        if not command:
            parser.error("run needs a command, after --")
        print()
        document["run"] = run_workload(readable, command, args.period)
        print_run(document["run"])
    if args.json:
        with open(args.json, "w") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
