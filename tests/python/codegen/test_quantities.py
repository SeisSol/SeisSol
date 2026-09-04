# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Unit tests for the quantity group descriptor.

The interesting cases are the four layouts SeisSol actually uses -- acoustic,
elastic, poroelastic and elastic-with-relaxation -- because the derived
sparsity and selectors have to reproduce what the equation classes previously
spelled out by hand.
"""

import numpy as np
import pytest
from kernels.quantities import (
    FaceRole,
    QuantityGroup,
    QuantityKind,
    layout,
    role_offset,
    rotation_spp,
    total_extent,
    traction_selector,
    velocity_selector,
    well_formed,
)

ACOUSTIC = [
    QuantityGroup("pprime", QuantityKind.SCALAR, FaceRole.TRACTION),
    QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
]

ELASTIC = [
    QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
    QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
]

POROELASTIC = ELASTIC + [
    QuantityGroup("p", QuantityKind.SCALAR, FaceRole.EXTRA_TRACTION),
    QuantityGroup("vf", QuantityKind.VECTOR, FaceRole.EXTRA_VELOCITY),
]

ELASTIC_MECHANISM = [QuantityGroup("theta", QuantityKind.SYM_TENSOR2)]
ACOUSTIC_MECHANISM = [QuantityGroup("theta", QuantityKind.SCALAR)]


class TestLayout:
    def test_extents(self):
        assert total_extent(layout(ACOUSTIC)) == 4
        assert total_extent(layout(ELASTIC)) == 9
        assert total_extent(layout(POROELASTIC)) == 13

    @pytest.mark.parametrize("mechanisms", [0, 1, 3, 5])
    def test_mechanism_repetitions(self, mechanisms):
        blocks = layout(ELASTIC, ELASTIC_MECHANISM, mechanisms)
        assert total_extent(blocks) == 9 + 6 * mechanisms
        blocks = layout(ACOUSTIC, ACOUSTIC_MECHANISM, mechanisms)
        assert total_extent(blocks) == 4 + mechanisms

    def test_offsets_are_contiguous(self):
        blocks = layout(POROELASTIC, ELASTIC_MECHANISM, 2)
        offset = 0
        for block in blocks:
            assert block.offset == offset
            offset += block.extent

    def test_mechanism_index_is_recorded(self):
        blocks = layout(ELASTIC, ELASTIC_MECHANISM, 3)
        assert [b.mechanism for b in blocks] == [None, None, 0, 1, 2]

    def test_negative_repetitions_rejected(self):
        with pytest.raises(ValueError):
            layout(ELASTIC, ELASTIC_MECHANISM, -1)


class TestRotationSparsity:
    def test_elastic_blocks(self):
        spp = rotation_spp(layout(ELASTIC))
        expected = np.zeros((9, 9), dtype=bool)
        expected[0:6, 0:6] = True
        expected[6:9, 6:9] = True
        assert np.array_equal(spp, expected)

    def test_acoustic_pressure_is_a_scalar(self):
        spp = rotation_spp(layout(ACOUSTIC))
        expected = np.zeros((4, 4), dtype=bool)
        expected[0, 0] = True
        expected[1:4, 1:4] = True
        assert np.array_equal(spp, expected)

    def test_poroelastic_four_blocks(self):
        spp = rotation_spp(layout(POROELASTIC))
        expected = np.zeros((13, 13), dtype=bool)
        expected[0:6, 0:6] = True
        expected[6:9, 6:9] = True
        expected[9, 9] = True
        expected[10:13, 10:13] = True
        assert np.array_equal(spp, expected)

    def test_viscoelastic_repeats_the_mechanism_block(self):
        spp = rotation_spp(layout(ELASTIC, ELASTIC_MECHANISM, 3))
        expected = np.zeros((27, 27), dtype=bool)
        expected[0:6, 0:6] = True
        expected[6:9, 6:9] = True
        for mech in range(3):
            offset = 9 + 6 * mech
            expected[offset : offset + 6, offset : offset + 6] = True
        assert np.array_equal(spp, expected)

    def test_viscoacoustic_memory_variables_are_scalars(self):
        spp = rotation_spp(layout(ACOUSTIC, ACOUSTIC_MECHANISM, 3))
        # Seven quantities: one pressure, three velocities, three scalar
        # memory variables -- every one of the latter rotates on its own.
        assert spp.shape == (7, 7)
        assert np.array_equal(np.diag(spp), np.ones(7, dtype=bool))
        assert spp[1:4, 1:4].all()
        assert not spp[4:7, 0:4].any()
        assert spp[4:7, 4:7].sum() == 3

    def test_blocks_are_disjoint_and_cover_everything(self):
        for groups, mech, reps in [
            (ACOUSTIC, ACOUSTIC_MECHANISM, 2),
            (ELASTIC, ELASTIC_MECHANISM, 3),
            (POROELASTIC, (), 0),
        ]:
            blocks = layout(groups, mech, reps)
            spp = rotation_spp(blocks)
            assert spp.diagonal().all()
            assert spp.sum() == sum(block.extent**2 for block in blocks)


class TestSelectors:
    def test_elastic(self):
        blocks = layout(ELASTIC)
        expected_t = np.zeros((3, 9))
        expected_t[0, 0] = expected_t[1, 3] = expected_t[2, 5] = 1
        expected_v = np.zeros((3, 9))
        expected_v[0, 6] = expected_v[1, 7] = expected_v[2, 8] = 1
        assert np.array_equal(traction_selector(blocks), expected_t)
        assert np.array_equal(velocity_selector(blocks), expected_v)

    def test_acoustic_has_no_shear_traction(self):
        blocks = layout(ACOUSTIC)
        expected_t = np.zeros((3, 4))
        expected_t[0, 0] = 1
        expected_v = np.zeros((3, 4))
        expected_v[0, 1] = expected_v[1, 2] = expected_v[2, 3] = 1
        assert np.array_equal(traction_selector(blocks), expected_t)
        assert np.array_equal(velocity_selector(blocks), expected_v)

    def test_poroelastic_appends_one_row_each(self):
        blocks = layout(POROELASTIC)
        expected_t = np.zeros((4, 13))
        expected_t[0, 0] = expected_t[1, 3] = expected_t[2, 5] = expected_t[3, 9] = 1
        expected_v = np.zeros((4, 13))
        expected_v[0, 6] = expected_v[1, 7] = expected_v[2, 8] = expected_v[3, 10] = 1
        assert np.array_equal(traction_selector(blocks), expected_t)
        assert np.array_equal(velocity_selector(blocks), expected_v)

    def test_mechanisms_widen_but_do_not_populate(self):
        plain = layout(ELASTIC)
        with_mech = layout(ELASTIC, ELASTIC_MECHANISM, 3)
        for select in (traction_selector, velocity_selector):
            wide, narrow = select(with_mech), select(plain)
            assert wide.shape == (3, 27)
            assert np.array_equal(wide[:, :9], narrow)
            assert not wide[:, 9:].any()


class TestRoleOffset:
    @pytest.mark.parametrize(
        ("groups", "expected"), [(ACOUSTIC, 1), (ELASTIC, 6), (POROELASTIC, 6)]
    )
    def test_velocity_offset(self, groups, expected):
        assert role_offset(layout(groups), FaceRole.VELOCITY) == expected

    def test_missing_role_raises(self):
        blocks = layout([QuantityGroup("q", QuantityKind.SCALAR)])
        with pytest.raises(ValueError):
            role_offset(blocks, FaceRole.VELOCITY)


class TestWellFormed:
    def test_accepts_the_real_layouts(self):
        assert well_formed(layout(ACOUSTIC), 4)
        assert well_formed(layout(ELASTIC), 9)
        assert well_formed(layout(POROELASTIC), 13)
        assert well_formed(layout(ELASTIC, ELASTIC_MECHANISM, 3), 27)

    def test_rejects_a_wrong_quantity_count(self):
        with pytest.raises(ValueError, match="covers 9 quantities"):
            well_formed(layout(ELASTIC), 13)

    def test_rejects_duplicate_names(self):
        groups = [
            QuantityGroup("q", QuantityKind.SCALAR),
            QuantityGroup("q", QuantityKind.VECTOR),
        ]
        with pytest.raises(ValueError, match="duplicate"):
            well_formed(layout(groups))

    def test_rejects_two_velocity_groups(self):
        groups = ELASTIC + [QuantityGroup("vf", QuantityKind.VECTOR, FaceRole.VELOCITY)]
        with pytest.raises(ValueError, match="at most one velocity"):
            well_formed(layout(groups))
