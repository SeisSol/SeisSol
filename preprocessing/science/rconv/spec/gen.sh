#!/bin/bash
ncgen -l c nrf.cdl > nrf.c
sed -i 's/Vector3_PERIOD_//g' nrf.c
sed -i 's/Subfault_units_PERIOD_//g' nrf.c
sed -i 's/Subfault_PERIOD_//g' nrf.c
