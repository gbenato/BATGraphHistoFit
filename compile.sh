g++ -std=c++11 `root-config --cflags --glibs` macros/BAT_Efficiency.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/BAT_Efficiency `bat-config --libs --cflags`
g++ -std=c++11  macros/PrintEfficiencyValues.cxx -o bin/PrintEfficiencyValues `root-config --cflags --glibs`
g++ -std=c++11  macros/CombineEfficiency.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/CombineEfficiency `bat-config --libs --cflags`
g++ -std=c++11 `root-config --cflags --glibs` macros/BAT_Reso.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/BAT_Reso `bat-config --libs --cflags`
g++ -std=c++11  macros/CombineReso.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/CombineReso `bat-config --libs --cflags`
