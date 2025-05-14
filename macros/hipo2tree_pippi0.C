#include "../src/CutManager.C"
#include "../src/CLAS12Ana.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"

//clas12root -l -b -q 'macros/hipo2tree_pippi0.C("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/nSidis_005032.hipo" , "out/test_pippi0/nSidis_005032.root", 500)'

int hipo2tree_pippi0(const char* hipoFile = "",
    const char* outputFile = "",
    const int maxEvents = 100){

    // Create a TFile to save the data
    TFile* fOut = new TFile(outputFile, "RECREATE");
    // Create a TTree to store the data
    EventTree* treeReco = new EventTree("EventTree");
    
    // Configure CLAS12 Reader and HipoChain
    // -------------------------------------
    clas12root::HipoChain _chain;
    _chain.Add(hipoFile);
    //_chain.Add("hipodir/nSidis_*") //for adding more files

    //use the config to select only the events in the hipo files that have the particles we are interested in
    auto _config_c12=_chain.GetC12Reader();
    _config_c12->addAtLeastPid(11, 1); //make sure there is at least one electron
    _config_c12->addAtLeastPid(211, 1); // at least one pi plus
    _config_c12->addAtLeastPid(22, 2); // at least 2 photons (pi0 decay products)
    //_config_c12->addAtLeastPid(2112,1); //at least one DETECTED neutron

    //generate the object that contains all the information for the desired events
    auto& _c12 = _chain.C12ref();
    
    // Check if "rg-a" is present in the hipoFile path
    bool doQADB = false;
    std::string hipoFilePath = hipoFile;
    if (hipoFilePath.find("rg-a") != std::string::npos &&
        hipoFilePath.find("montecarlo") == std::string::npos) {
        doQADB = true;
    }
    //run a quality assurance script made for this run data set
    if (doQADB) {
        _config_c12->applyQA("pass1");
        _config_c12->db()->qadb_addQARequirement("OkForAsymmetry");
    }

    fOut->Close();
    return 0;
}