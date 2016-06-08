# Move to the project root.
BASEDIR=$(dirname "$0")
cd $BASEDIR

# Create the debug folder.
rm -rf ./debug
mkdir -p ./debug

# Compile the CUDA program.
nvcc -c --std=c++11 -arch=sm_30 ./src/time_series_aaft.cu -o ./debug/time_series_aaft.o
nvcc -c --std=c++11 -arch=sm_30 ./src/fmri_corr_coef.cu -o ./debug/fmri_corr_coef.o
nvcc --std=c++11 -arch=sm_30 ./debug/time_series_aaft.o ./debug/fmri_corr_coef.o ./src/main.cu -o ./debug/out

# Execute.
cd ./debug


# AD1 Left
./out "../lol_stc_csv/s01_120415-002-AD1-lh.csv" "../lol_stc_csv/s02_120415-002-AD1-lh.csv" "../lol_stc_csv/s03_012116-004-AD1-lh.csv" "../lol_stc_csv/s09_012816-004-AD1-lh.csv" "../lol_stc_csv/s11_032516-003-AD1-lh.csv" "../lol_stc_csv/s14_032816-002-AD1-lh.csv" "../lol_stc_csv/s17_042416-002-AD1-lh.csv"
# AD1 Right
# ./out "../lol_stc_csv/s01_120415-002-AD1-rh.csv" "../lol_stc_csv/s02_120415-002-AD1-rh.csv" "../lol_stc_csv/s03_012116-004-AD1-rh.csv" "../lol_stc_csv/s09_012816-004-AD1-rh.csv" "../lol_stc_csv/s11_032516-003-AD1-rh.csv" "../lol_stc_csv/s14_032816-002-AD1-rh.csv" "../lol_stc_csv/s17_042416-002-AD1-rh.csv"


# SUP1 Left
# ./out "../lol_stc_csv/s11_032516-002-SU1-lh.csv" "../lol_stc_csv/s12_032516-004-SU1-lh.csv" "../lol_stc_csv/s13_032816-002-SU1-lh.csv" "../lol_stc_csv/s14_032816-003-SU1-lh.csv" "../lol_stc_csv/s15_042416-009-SU1-lh.csv" "../lol_stc_csv/s16_042416-004-SU1-lh.csv" "../lol_stc_csv/s17_042416-005-SU1-lh.csv"
# SUP1 Right
# ./out "../lol_stc_csv/s11_032516-002-SU1-rh.csv" "../lol_stc_csv/s12_032516-004-SU1-rh.csv" "../lol_stc_csv/s13_032816-002-SU1-rh.csv" "../lol_stc_csv/s14_032816-003-SU1-rh.csv" "../lol_stc_csv/s15_042416-009-SU1-rh.csv" "../lol_stc_csv/s16_042416-004-SU1-rh.csv" "../lol_stc_csv/s17_042416-005-SU1-rh.csv"


# AD2 Left
# ./out "../lol_stc_csv/s01_120415-006-AD2-lh.csv" "../lol_stc_csv/s02_120415-006-AD2-lh.csv" "../lol_stc_csv/s03_012116-003-AD2-lh.csv" "../lol_stc_csv/s04_012116-003-AD2-lh.csv" "../lol_stc_csv/s05_012116-003-AD2-lh.csv" "../lol_stc_csv/s06_012116-003-AD2-lh.csv" "../lol_stc_csv/s07_012116-004-AD2-lh.csv" "../lol_stc_csv/s08_012816-004-AD2-lh.csv" "../lol_stc_csv/s09_012816-003-AD2-lh.csv" "../lol_stc_csv/s10_012816-004-AD2-lh.csv"
# AD2 Right
# ./out "../lol_stc_csv/s01_120415-006-AD2-rh.csv" "../lol_stc_csv/s02_120415-006-AD2-rh.csv" "../lol_stc_csv/s03_012116-003-AD2-rh.csv" "../lol_stc_csv/s04_012116-003-AD2-rh.csv" "../lol_stc_csv/s05_012116-003-AD2-rh.csv" "../lol_stc_csv/s06_012116-003-AD2-rh.csv" "../lol_stc_csv/s07_012116-004-AD2-rh.csv" "../lol_stc_csv/s08_012816-004-AD2-rh.csv" "../lol_stc_csv/s09_012816-003-AD2-rh.csv" "../lol_stc_csv/s10_012816-004-AD2-lh.csv"


# SUP2 Left
# ./out "../lol_stc_csv/s04_012116-004-SU2-lh.csv" "../lol_stc_csv/s05_012116-004-SU2-lh.csv" "../lol_stc_csv/s06_012116-004-SU2-lh.csv" "../lol_stc_csv/s07_012116-003-SU2-lh.csv" "../lol_stc_csv/s08_012816-003-SU2-lh.csv" "../lol_stc_csv/s10_012816-003-SU2-lh.csv" "../lol_stc_csv/s12_032516-005-SU2-lh.csv" "../lol_stc_csv/s13_032816-003-SU2-lh.csv" "../lol_stc_csv/s15_042416-007-SU2-lh.csv"
# SUP2 Right
# ./out "../lol_stc_csv/s04_012116-004-SU2-rh.csv" "../lol_stc_csv/s05_012116-004-SU2-rh.csv" "../lol_stc_csv/s06_012116-004-SU2-rh.csv" "../lol_stc_csv/s07_012116-003-SU2-rh.csv" "../lol_stc_csv/s08_012816-003-SU2-rh.csv" "../lol_stc_csv/s10_012816-003-SU2-rh.csv" "../lol_stc_csv/s12_032516-005-SU2-rh.csv" "../lol_stc_csv/s13_032816-003-SU2-rh.csv" "../lol_stc_csv/s15_042416-007-SU2-rh.csv"

