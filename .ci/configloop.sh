#! /usr/bin/env sh
# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

configloop() {
    local function="$1"
    local cfgfile="$2"

    SUMMARY="Summary:"
    FAILURE=0

    COUNT=$( jq '. | length' "$cfgfile" )
    echo "Iterating over $COUNT configs..."
    for (( i=0; i<$COUNT; i++ )) ; do
        CONFIG=$( jq -rc ".[$i]" "$cfgfile" )
        echo "--------------------------------------------------------------------------------"
        echo "Current config ($i): $CONFIG"

        mkdir -p build-$i
        cd build-$i

        set +e

        $function $i $CONFIG

        RETVAL=$?

        set -e

        echo "Result ($i): $RETVAL"
        echo "--------------------------------------------------------------------------------"

        SUMMARY="$SUMMARY\n$CONFIG ($i): $RETVAL"

        if [ $RETVAL != 0 ]; then
            FAILURE=1
        fi

        cd ..
    done

    echo $SUMMARY

    exit $FAILURE
}

testloop() {
    local function="$1"
    local cfgfile="$2"

    SUMMARY="Summary:"
    FAILURE=0

    COUNT=$( jq '. | length' "$cfgfile" )
    echo "Iterating over $COUNT configs..."
    for (( i=0; i<$COUNT; i++ )) ; do
        CONFIG=$( jq -rc ".[$i]" "$cfgfile" )
        echo "--------------------------------------------------------------------------------"
        echo "Current config ($i): $CONFIG"

        set +e

        $function $i $CONFIG

        RETVAL=$?

        set -e

        echo "Result ($i): $RETVAL"
        echo "--------------------------------------------------------------------------------"

        SUMMARY="$SUMMARY\n$CONFIG ($i): $RETVAL"

        if [ $RETVAL != 0 ]; then
            FAILURE=1
        fi
    done

    echo $SUMMARY

    exit $FAILURE
}

configname() {
    local CONFIG="$1"
    local DEVICE="$2"
    ORDER=$(echo "$CONFIG" | jq -rc '.order')
    EQUATION=$(echo "$CONFIG" | jq -rc '.equation')
    MECHANISMS=$(echo "$CONFIG" | jq -rc '.mechanisms')
    PRECISION=$(echo "$CONFIG" | jq -rc '.precision')
    SIMULATIONS=$(echo "$CONFIG" | jq -rc '.simulations')

    if [ $PRECISION = single ]; then
        PRECISIONPFX=s;
    else
        PRECISIONPFX=d;
    fi;
    if [ $SIMULATIONS == "1" ]; then
        MULTISIM='';
    else
        MULTISIM="-f$SIMULATIONS";
    fi;

    EXENAME="$DEVICE-$EQUATION-$PRECISIONPFX-p$ORDER$MULTISIM"

    echo $EXENAME
}
