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
        _config_c12->applyQA("pass1"); //specifies which QADB to use
        //_config_c12->db()->qadb_addQARequirement("OkForAsymmetry"); //specifies which events to cut from the data. It seems that "OkForAsymmetry" was discontinued: https://github.com/JeffersonLab/clas12-qadb?tab=readme-ov-file#info
    }

    // Add Analysis Objects
    // -------------------------------------
    CutManager _cm = CutManager();
    CLAS12Ana clas12ana = CLAS12Ana(_c12);
    clas12ana.set_run_config(_c12); //grabs the event number index, run number index, and torus info index, saving it to the corresponding object private methods for clas12ana

    // Add Analysis Structs - from the Structs.h src file
    // -------------------------------------
    std::vector<part> vec_particles;
    std::vector<part> vec_mcparticles;
    EVENT_INFO event_info;
    EVENT event;

    int whileidx = 0;
    int _ievent = 0;
    int badAsym = 0;

    //now loop over the chain of events (chain becomes like on big hipofile) until we reach the max number of events. note: if maxEvents = -1, then it will run through the whole chain.
    while (_chain.Next() == true && (whileidx < maxEvents || maxEvents < 0)) {
        //this if statement just prints to the console every 10000 events. It also mentions if any events are skipped from the quality assurance check from QADB.
        if (whileidx % 10000 == 0 && whileidx != 0) {
            std::cout << whileidx << " events read | "
                      << _ievent * 100.0 / whileidx << "% passed event selection | "
                      << badAsym << " events skipped from QADB"
                      << std::endl;
        }

        //for each event, populate the neccesary information in the clas12ana and cutmanager contexts using info about each event found in the constants class
        //----------------------------------------------------------------------------------------------------------------------------------------------------------
        clas12ana.get_event_info(_c12, event_info); //this line fills event_info with the relevant properties of that event based found in constants.h
        event_info.uID = whileidx; //keeps track of the iteration in the event_info method
        whileidx++;
        _cm.set_run(event_info.run); //sets relevant runinfo for the cutmanager class to use


        //Now perform cuts
        //----------------------------------------------------------------------------------------------------------------------------------------------------------
        vec_particles = clas12ana.load_reco_particles(_c12); //populates the initial vector of particles to be cut from
        //determine which particle is the scattered electron
        int idx_scattered_ele = clas12ana.find_reco_scattered_electron(vec_particles);
        if (idx_scattered_ele == -1)
            continue; // No scattered electron found in FD
        vec_particles[idx_scattered_ele].is_scattered_electron = 1;

        clas12ana.fill_reco_event_variables(event, vec_particles);
        if (event.y > 0.8 || event.Q2 < 1||event.W<2) //DIS cuts
            continue; 
        vec_particles = _cm.filter_particles(vec_particles); // Apply all other Cuts

        //final check that this event is meaningful
        int num_e = 0;
        int num_piplus = 0;
        int num_piminus = 0;
        int num_gamma = 0;
        int num_proton = 0;
        for (auto part : vec_particles) {
            if (part.pid == 11) num_e++;
            if (part.pid == 211) num_piplus++;
            if (part.pid == -211) num_piminus++;
            if (part.pid == 22) num_gamma++;
            if (part.pid == 2212) num_proton++;
        }
        if (num_e == 0 || num_piplus < 1 || num_gamma < 2) continue;

        //then fill the tree with the vector of particles, event, and event_info
        treeReco->FillTree(vec_particles, event, event_info);
        _ievent++;
    }

    fOut->cd();
    treeReco->Write();
    fOut->Close();
    return 0;
}