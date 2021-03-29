if [ $# -ne 1 ]; then
    echo "$0: Takes one argument, the name of the AI to compile"
    exit 1
fi
AI_Name=$1 #Get the AI name you passed as a command line parameter
AI_Name="${AI_Name%.*}" #Remove extension in case you passed "V4.cpp" instead of "V4"

clang++-9 -std=c++17 -march=native -mpopcnt -mbmi2 -mfma -mavx2 -Ofast -funroll-loops -finline "$AI_Name.cpp" -lpthread -o "$AI_Name"

