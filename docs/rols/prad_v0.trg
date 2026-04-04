 
arch = Red

Following environment variables are defined:
OSTYPE          = Linux
PLATFORM        = x86_64
MACHINE         = x86_64_RHEL9
OSTYPE_PLATFORM = Linux_x86_64
OSTYPE_MACHINE  = Linux_x86_64_RHEL9

#
# main config file
#
#include ti/clasdev.cnf
#include dsc2/tdcpcal2.cnf
#include fadc250/adcpcal2.cnf
#include tdc1190/clasdev.cnf

include fadc250/peds/adchycal1_ped.cnf
include fadc250/peds/adchycal2_ped.cnf
include fadc250/peds/adchycal3_ped.cnf
include fadc250/peds/adchycal4_ped.cnf
include fadc250/peds/adchycal5_ped.cnf
include fadc250/peds/adchycal6_ped.cnf
include fadc250/peds/adchycal7_ped.cnf


### Align timing of crates
TI_CRATE adchycal1
TI_FIBER_DELAY_OFFSET 0x80 0xcf
TI_CRATE end
TI_CRATE adchycal2
TI_FIBER_DELAY_OFFSET 0x80 0xcc
TI_CRATE end
TI_CRATE adchycal3
TI_FIBER_DELAY_OFFSET 0x80 0xd3
TI_CRATE end
TI_CRATE adchycal4
TI_FIBER_DELAY_OFFSET 0x80 0xd0
TI_CRATE end
TI_CRATE adchycal5
TI_FIBER_DELAY_OFFSET 0x80 0xcc
TI_CRATE end
TI_CRATE adchycal6
TI_FIBER_DELAY_OFFSET 0x80 0xd0
TI_CRATE end
TI_CRATE adchycal7
TI_FIBER_DELAY_OFFSET 0x80 0xd0
TI_CRATE end

VTP_CRATE adchycal1vtp
#      fiber:  1  2  3  4
VTP_FIBER_EN   1  1  0  0
VTP_CRATE end

VTP_CRATE adchycal2vtp
VTP_FIBER_EN   1  1  1  1
VTP_CRATE end

VTP_CRATE adchycal3vtp
VTP_FIBER_EN   1  1  1  1
VTP_CRATE end

VTP_CRATE adchycal4vtp
VTP_FIBER_EN   1  1  1  1
VTP_CRATE end

VTP_CRATE adchycal5vtp
VTP_FIBER_EN   1  1  1  1
VTP_CRATE end

VTP_CRATE adchycal6vtp
VTP_FIBER_EN   1  1  1  1
VTP_CRATE end

VTP_CRATE adchycal7vtp
VTP_FIBER_EN   1  1  0  0
VTP_CRATE end


TS_HOLDOFF   1  5 1
TS_HOLDOFF   2  5 1
TS_HOLDOFF   3 15 1
TS_HOLDOFF   4 10 1


TS_CRATE trig0

TS_BLOCK_LEVEL   1
TS_BUFFER_LEVEL  1
#TS_BLOCK_LEVEL   5
#TS_BUFFER_LEVEL  5


#  FP MASK
#  0x0000FF00 - SSP
#  0x00000100 - PRAD TRGBIT 0
#  0x00000200 - PRAD TRGBIT 1
#  0x00000400 - PRAD TRGBIT 2
#  0x00000800 - PRAD TRGBIT 3
#  0x00001000 - PRAD TRGBIT 4
#  0x00002000 - PRAD TRGBIT 5
#  0x00004000 - PRAD TRGBIT 6
#  0x00008000 - PRAD TRGBIT 7

#   V1495
# 0xFFFF0000 - externals
# 0x00800000 - OR from SD/FADC  
# 0x01000000 - LMS
# 0x02000000 - alpha
# 0x04000000 - Faraday
# 0x08000000 - Master OR


#TS_FP_INPUT_MASK  0x00000100  # SSP PRAD TRGBIT0
#TS_FP_INPUT_MASK  0x80000100  # SSP PRAD TRGBIT0 OR trig15 from v1495
#TS_FP_INPUT_MASK  0x80000000  # trig15 from v1495

TS_FP_INPUT_MASK  0x0F00FF00  # SSP PRAD TRGBIT 0-7 LMS alpha Faraday Master-OR
#TS_FP_INPUT_MASK  0x0F000000  # 

###    Prescale (1+2^(x-1) for x>0; x=0 prescale = 1, no prescale
###    Prescales    0:1 1:2   2:3   3:5   4:9  5:17 15:16385

TS_FP_PRESCALE  24     0   # FP TriggerBit24: LMS
TS_FP_PRESCALE  25    15   # FP TriggerBit25: alpha
TS_FP_PRESCALE  26    15   # FP TriggerBit26: Faraday
TS_FP_PRESCALE  27    15   # FP TriggerBit27: Master-OR

#   Random trigger run ONLY
#   0 - disable
#   1 - enable; 
#   Second number is frequency: =500KHz/2^x  4:31KHz  5:15KHz 6:7.8KHz 7:3.9KHz  15:15 Hz 
#TS_RANDOM_TRIGGER 1   5
TS_RANDOM_TRIGGER 0   5


TS_CRATE end


###############################################
SSP_CRATE all
SSP_SLOT all
SSP_W_OFFSET   2700
SSP_W_WIDTH    400

SSP_SET_IO_SRC  16   5   #SYNC=FPLVDSIN0
SSP_SET_IO_SRC  15   6   #TRIG=FPLVDSIN1
SSP_SET_IO_SRC   7   20  #P2OUT0=TRGBIT0
SSP_SET_IO_SRC   8   21  #P2OUT1=TRGBIT1
SSP_SET_IO_SRC   9   22  #P2OUT2=TRGBIT2
SSP_SET_IO_SRC  10   23  #P2OUT3=TRGBIT3
SSP_SET_IO_SRC  11   24  #P2OUT4=TRGBIT4
SSP_SET_IO_SRC  12   25  #P2OUT5=TRGBIT5
SSP_SET_IO_SRC  13   26  #P2OUT6=TRGBIT6
SSP_SET_IO_SRC  14   18  #P2OUT7=100Hz Pulser TRGBIT7

# SSP pulser rate in Hz
SSP_PULSER  100

#                   Latency(ns)
#                   |    TrgWidth(ns)
#                   |    |
SSP_PRAD_TRIGGER    1540 100

#               N[0-7]
#               |    Prescale
#               |    |    ClusterMultDelay(ns)
#               |    |    |    ClusterSumDelay(ns)
#               |    |    |    |    RawSumDelay(ns)
#               |    |    |    |    |    ClusterMultMin
#               |    |    |    |    |    |    ClusterSumMin(MeV-when calibrated w/gain)
#               |    |    |    |    |    |    |     RawSumMin
#               |    |    |    |    |    |    |     | Flags0
#               |    |    |    |    |    |    |     | | Flags1
#               |    |    |    |    |    |    |     | | | ShortName
#               |    |    |    |    |    |    |     | | | |         LongName
#               |    |    |    |    |    |    |     | | | |         |
SSP_PRAD_TRGBIT 0    1    0    0    0    0    0  1000 0 0 Trgbit0   Trgbit0 # 0: Raw sum min > 1000 
SSP_PRAD_TRGBIT 1    0    0    0    0    1 1000     0 0 0 Trgbit1   Trgbit1 # 1: >=1 cluster with cluster sum>1000MeV
SSP_PRAD_TRGBIT 2    0    0    0    0    2 1000     0 0 0 Trgbit2   Trgbit2 # 2: >=2 cluster with cluster sum>1000MeV
SSP_PRAD_TRGBIT 3    0    0    0    0    3 1000     0 0 0 Trgbit3   Trgbit3 # 3, >=3 cluster with cluster sum>1000MeV
SSP_PRAD_TRGBIT 4    0    0    0    0  255    1     0 0 0 Trgbit4   Trgbit4 # 4, 
SSP_PRAD_TRGBIT 5    0    0    0    0  255    1     0 0 0 Trgbit5   Trgbit5 # 5, 
SSP_PRAD_TRGBIT 6    0    0    0    0  255    1     0 0 0 Trgbit6   Trgbit6 # 6, 


SSP_CRATE end
###############################################


###############################################
VTP_CRATE all
VTP_W_OFFSET 3000
VTP_W_WIDTH  400

VTP_PRAD_FADCSUM_CH 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF 0xFFFF

#        slot: 10 13  9 14  8 15  7 16  6 17  5 18  4 19  3 20
#     payload:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
VTP_PAYLOAD_EN  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1

#                 SeedThr(MeV)
#                 |    ClusterHitDt(+/-ns, max +/-16ns)
#                 |    |    ClusterNHitsMin
#                 |    |    |    ClusterThr(MeV)
#                 |    |    |    |    EnergySumIntWidth(ns)
#                 |    |    |    |    |
VTP_PRAD_CLUSTER 100  16    2  200   32     # 1MeV seed thr, +/-16ns cluster hit coincidence, 1 hit min in cluster, 1MeV cluster threshold, 32ns integration width for raw and cluster energy sum

VTP_CRATE end


