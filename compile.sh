g++ -std=c++17 `root-config --cflags --glibs` macros/BAT_Efficiency.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/BAT_Efficiency `bat-config --libs --cflags`
g++ -std=c++17 `root-config --cflags --glibs` macros/CombineEfficiency.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/CombineEfficiency `bat-config --libs --cflags`
g++ -std=c++17 `root-config --cflags --glibs` macros/BAT_Reso.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/BAT_Reso `bat-config --libs --cflags`
