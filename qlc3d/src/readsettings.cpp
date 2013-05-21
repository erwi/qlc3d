#include <stdio.h>
#include <stdlib.h>
#include <qlc3d.h>
#include <string>
#include <vector>
#include <reader.h>
#include <iostream>
#include <meshrefinement.h>
#include <refinfo.h>
#include <filesysfun.h>
using namespace std;


void problem(const char *error)
{
    /*! Prints warning message for a bad setting and exits*/
    printf("WARNING - problem reading setting for : ");
    printf("%s\n",error);
    printf("\n");
    exit(1);
}
void problem(std::string& name, int ret){
    /*! Checks reader return value. If it is an error, prints error message and quits.
        Use this for mandatory settings.
     */

    if ( ret == READER_SUCCESS )
        return;

    Reader temp;
    if (ret!= READER_SUCCESS){
        cout << "Problem reading " << name << " , computer says: " << temp.getErrorString(ret) << endl;
        cout << "bye!"<<endl;
        exit(1);
    }
}
void problemo(std::string& name, int ret)
{
    /*! Checks reader return value. If it is an error, prints warning message but does not quit.
       Use this for optional settings.*/
    Reader temp;
    if ( (ret!= READER_SUCCESS) )
    {
        cout << "Problem reading " << name << " , computer says: " << temp.getErrorString(ret) << endl;
    }
}
// end problemo (optional parameter problem)
void problem_format(std::string &name, int ret)
{
    // CHECKS READER RETURN VALUE. QUITS IN CASE OF BAD FORMAT, BUT NOT IF
    // VALUE WAS NOT FOUND

    // RETURN IF ...
    if  ( (ret == READER_NOT_FOUND) ||  // NOT FOUND
          (ret == READER_SUCCESS) )     // EVERYTHING WENT BETTER THAN EXPECTED
        return;

    Reader temp; // SHOULD MAKE A STATIC FUNCTION INSTEAD
    std::string error = temp.getErrorString( ret );
    cout << "Problem reading " << name << " , computer says: " << error << endl;
    cout << "bye!"<<endl;
    exit(1);
}

bool isOK(const std::string& name, int ret)
{
    // EXITS AND PRINTS ERROR IF PROBLEM READING A VALUE THAT IS FOUND.
    // USE THIS TO CHECK THE SYNTAX OF OPTIONAL VALUE DEFINITIONS

    // ALL OK.
    if ( ret == READER_SUCCESS ) return true;


    if ( ( ret == READER_BAD_FORMAT ) ||
         (ret == READER_BAD_VALUE  ) ||
         (ret == READER_BAD_FILE) )
    {
        Reader temp;
        cout << "Problem reading " << name << ", computer says: " << temp.getErrorString( ret ) << endl;
        exit(ret);
    }

    return false;   // SILENTLY RETURN FALSE IF KEY/VALUE PAIR WAS NOT FOUND

}


std::string setStructureKey(const char* struct_name, const unsigned int number, const char* key_name)
{
    std::stringstream ss;
    ss << struct_name << number <<"."<<key_name;
    std::string key;
    ss >> key;
    return key;
}


void readLC(LC& lc,Reader& reader)
{

    /// THIS CHECKING SHOULD BE DONE IN THE READER CLASS
    //	reader.file.seekg(0);
    //	reader.file.clear();
    if ( ! reader.file.good() ){
        cout << "error reading LC - bye!" <<endl;
        exit(1);
    }
    string name = "";
    int ret     = 0;
    double val  = 0;
    // ELASTIC COEFFS
    name = "K11";
    ret = reader.readNumber(name , val);
    if (ret == READER_SUCCESS)
        lc.K11 = val;
    problem(name, ret);

    name = "K22";
    ret = reader.readNumber(name , val);
    if (ret == READER_SUCCESS)
        lc.K22 = val;
    problem(name, ret);

    name = "K33";
    ret = reader.readNumber(name , val);
    if (ret == READER_SUCCESS)
        lc.K33 = val;
    problem(name, ret);

    name = "p0";
    ret = reader.readNumber(name , val);
    if(ret == READER_SUCCESS)
        lc.p0 = val;
    problem_format(name, ret);
    //problemo(name, ret);
    // THERMOTROPIC COEFFS
    name = "A";
    ret = reader.readNumber(name , val);
    if(ret == READER_SUCCESS)
        lc.A = val;
    problem(name,ret);

    name = "B";
    ret = reader.readNumber(name, val);
    if(ret == READER_SUCCESS)
        lc.B = val;
    problem(name,ret);

    name = "C";
    ret = reader.readNumber(name, val);
    if(ret == READER_SUCCESS)
        lc.C = val;
    problem(name,ret);

    // ELECTRIC COEFFS
    name = "eps_par";
    ret = reader.readNumber(name, val);
    if(ret == READER_SUCCESS)
        lc.eps_par = val;
    problem(name, ret);

    name = "eps_per";
    ret = reader.readNumber(name, val);
    if(ret == READER_SUCCESS)
        lc.eps_per = val;
    problem(name, ret);

    name = "e11";
    ret = reader.readNumber(name, val);
    problem_format(name, ret);
    if(ret == READER_SUCCESS)
        lc.e11 = val;


    name = "e33";
    ret = reader.readNumber(name, val);
    problem_format(name, ret);
    if(ret == READER_SUCCESS)
        lc.e33 = val;


    // VISCOUS COEFFICIENTS
    name = "gamma1";
    ret = reader.readNumber(name, val);
    if(ret == READER_SUCCESS)
        lc.gamma1 = val;
    problem(name, ret);

    name = "gamma2";
    ret = reader.readNumber(name, val);
    problem_format(name, ret);
    if(ret == READER_SUCCESS)
        lc.gamma2 = val;

}//end void readLC


void readSimu(Simu* simu, Reader& reader, EventList& evel)
{
    std::string name;
    std::string str_var;
    int         int_var = 0;
    double      dbl_var = 0;
    int         ret;

    //  READ STRING SETTINGS
    //  MANDATORY STRING SETTINGS

    ret = 0;
    name = "EndCriterion";
    ret = reader.readString(name , str_var);
    

    if ( ret== READER_SUCCESS)
        simu->setEndCriterion(str_var);

    problem(name, ret);

    name = "MeshName";
    ret = reader.readString(name, str_var);
    if ( ret == READER_SUCCESS )
        simu->MeshName = str_var;

    problem(name, ret);

    // OPTIONAL STRING SETTINGS
    name = "LoadQ";
    ret = reader.readString(name , str_var);
    if (ret == READER_SUCCESS)
        simu->setLoadQ(str_var);
    //problemo(name, ret);

    name = "SaveDir";
    ret = reader.readString(name , str_var);
    if ( ret == READER_SUCCESS) // if specified, use it
    {
        simu->setSaveDir( simu->getCurrentDir() + "/" + str_var);
    }
    else // otherwise use default ( = res )
    {
        simu->setSaveDir( simu->getCurrentDir() + "/" + simu->getSaveDir() );
    }

    name = "QMatrixSolver";
    ret = reader.readString(name, str_var);
    if (ret == READER_SUCCESS)
        simu->setQMatrixSolver(str_var);

    //========================
    // SCALAR VALUES
    //========================
    name = "SaveIter";
    ret = reader.readNumber(name, int_var);
    if (ret == READER_SUCCESS){
        simu->setSaveIter(int_var);
        evel.setSaveIter( (size_t) int_var );
    }

    name = "SaveTime";
    ret = reader.readNumber(name, dbl_var );
    if ( isOK(name, ret) ){
        evel.setSaveTime( dbl_var );
    }


    name = "EndValue";
    ret = reader.readNumber(name, dbl_var);
    if (ret == READER_SUCCESS){
        simu->setEndValue(dbl_var);
    }
    problem(name, ret);

    name = "dt";
    ret = reader.readNumber(name , dbl_var);
    if (ret == READER_SUCCESS){
        simu->setdt(dbl_var);
    }


    name = "TargetdQ";
    ret = reader.readNumber(name , dbl_var);
    if (ret == READER_SUCCESS){
        simu->setTargetdQ(dbl_var);
    }
    else{
        exit(1);
    }
    //problemo(name, ret);

    name = "Maxdt";
    ret = reader.readNumber(name , dbl_var);
    if(ret==READER_SUCCESS)
        simu->setMaxdt(dbl_var);
    //problemo(name, ret);

    name = "MaxError";
    ret = reader.readNumber(name, dbl_var);
    if(ret==READER_SUCCESS)
        simu->setMaxError(dbl_var);
    //problemo(name, ret);

    name = "OutputEnergy";
    ret  = reader.readNumber(name , int_var);
    if ( ret == READER_SUCCESS)
        simu->setOutputEnergy(int_var);
    //problemo(name, ret);

    name = "OutputFormat";
    ret = reader.readNumber(name , int_var);
    if(ret == READER_SUCCESS)
        simu->setOutputFormat(int_var);

    name = "NumAssemblyThreads";
    ret = reader.readNumber(name, int_var);
    if (ret == READER_SUCCESS)
        simu->setAsseblyThreadCount( (unsigned int) int_var);

    name = "NumMatrixSolverThreads";
    ret = reader.readNumber(name, int_var);
    if (ret == READER_SUCCESS)
        simu->setMatrixSolverThreadCount( (unsigned int) int_var);



    name = "SaveFormat";
    std::vector < std::string > vec_str;
    ret = reader.readStringArray( name , vec_str );

    if ( ret == READER_SUCCESS ){
        simu->clearSaveFormat(); // REMOVES DEFAULT
        for (size_t i = 0 ; i < vec_str.size() ; i++){
            simu->addSaveFormat( vec_str[i] );  // SAVE FORMATS ARE DEFINED IN Simu
        }
    }

    //------------------------------------------
    //
    // READ VECTOR VALUES
    //
    //------------------------------------------
    name = "StretchVector";
    vector <double> vec;

    ret = reader.readNumberArray(name , vec );
    if (ret == READER_SUCCESS){
        if (vec.size() != 3){
            cout << "error - invalid StretchVector length - bye!\n"<<endl;
            exit(1);
        }
        //cout << "stretch :"<< vec[0] <<","<< vec[1] <<"," << vec[2] << endl;
        simu->setStretchVectorX(vec[0]);
        simu->setStretchVectorY(vec[1]);
        simu->setStretchVectorZ(vec[2]);
        //cout << "stretch" << simu->getStretchVectorX() << "," << simu->getStretchVectorY() << "," << simu->getStretchVectorZ() << endl;
    }
    //problemo(name, ret);

    name = "EnergyRegion";
    ret = reader.readNumberArray(name , vec);
    if ((ret == READER_SUCCESS) && (vec.size() == 3)){
        simu->setEnergyRegionX(vec[0]);
        simu->setEnergyRegionY(vec[1]);
        simu->setEnergyRegionZ(vec[2]);
    }

    name = "dtLimits";
    ret = reader.readNumberArray(name , vec);
    if ( (ret == READER_SUCCESS) && (vec.size() == 2)){
        simu->setdtLimits(vec[0], vec[1]);
    }

    name = "dtFunction";
    ret = reader.readNumberArray(name , vec );
    if ((ret == READER_SUCCESS) && (vec.size() == 4) ){
        double F[4] ={vec[0] , vec[1] , vec[2], vec[3]};
        simu->setdtFunction(F);
    }

    name = "RegularGridSize";
    ret = reader.readNumberArray( name, vec);
    if (ret == READER_SUCCESS ){
        if (vec.size()!= 3){
            problem(name, READER_BAD_FORMAT );
        }
        simu->setRegularGridSize( (size_t) vec[0], (size_t) vec[1], (size_t) vec[2] );
    }


}//end void readSimu

void readBoxes(Boxes* boxes, Reader& reader) //config_t* cfg)
{
    int maxbox = 100; // maximum number of supported boxes...
    for (int i = 1 ; i < maxbox ; i++){ // for boxes 1-99
        // cout << i << endl;
        stringstream box;
        string name;
        string str_val;

        int ret;
        box << "BOX" << i <<".Type";
        box >> name;
        reader.file.seekg(0);
        reader.file.clear();
        if (!reader.file.good()){
            cout << "error reading Boxes, file pointer is bad - bye!" << endl;
            exit(1);
        }

        ret = reader.readString(name , str_val);
        //=============================
        // ADD BOXES IF THEY EXIST
        //=============================
        if (ret == READER_SUCCESS){ // if Boxi exists
            // cout << "reading box " << i << endl;
            Box* b = new Box(i);     // CREATE EMPTY BOX OBJECT
            b->setBoxType(str_val);
            //b->Type = str_val;

            vector <double> par;
            box.clear();
            name.clear();
            box << "BOX" << i <<".Params";
            box >> name;


            ret = reader.readNumberArray(name , par);
            problem(name, ret);
            b->setParams(par);

            // BOX X
            box.clear();
            name.clear();
            box << "BOX" << i <<".X";
            box >> name;
            ret = reader.readNumberArray(name , par);

            //cout << "X :" << par[0] << "," << par[1] << endl;
            problem(name, ret);
            b->setX(par);


            // BOX Y
            box.clear();
            name.clear();
            box << "BOX" << i <<".Y";
            box >> name;
            ret = reader.readNumberArray(name, par);
            problem(name, ret);
            b->setY(par);
            // BOX Z
            box.clear();
            name.clear();
            box << "BOX" << i <<".Z";
            box >> name;
            ret = reader.readNumberArray(name , par);
            problem(name , ret);
            b->setZ(par);
            // BOX TILT
            box.clear();
            name.clear();
            box << "BOX" << i <<".Tilt";
            box >> name;
            ret = reader.readNumberArray(name, par);
            problem(name , ret);
            b->setTilt(par);
            // BOX TWIST
            box.clear();
            name.clear();
            box << "BOX" << i <<".Twist";
            box >> name;
            ret = reader.readNumberArray(name, par);
            problem(name, ret);
            b->setTwist(par);

            // ADD TO LIST OF BOXES
            boxes->addBox(b);

        }// end if box exists
    }
    //end for i = 1:maxbox

}// end void readBoxes

void readAlignment(Alignment* alignment, Reader& reader)
{

    int num_surfaces = 99; // limit surfaces to 0-99
    std::string surface_name;
    std::string surface_setting_name;


    for ( int i = 0; i < num_surfaces; i++)	{
        stringstream ss;
        string name;
        string str_val;
        int ret = 0;
        ss << "FIXLC" << i << ".Anchoring";
        ss >> name;

        ret = reader.readString(name , str_val);

        if ( ret == READER_SUCCESS){
            Surface* s = new Surface(i);
            s->setAnchoringType(str_val);

            double dbl_val;

            ss.clear(); name.clear();
            ss << "FIXLC" << i << ".Strength";
            ss >> name;
            ret = reader.readNumber(name , dbl_val);
            if ( ! s->getisFixed() ) {          // if not a fixed surface, set  mandatory strength
                problem(name , ret);
                s->setStrength(dbl_val);
            }
            else{
                if (ret==READER_SUCCESS)    //SET OPTIONAL
                    s->setStrength(dbl_val);
                else
                    s->setStrength(99e99);
            }


            ss.clear(); name.clear();
            ss << "FIXLC"<<i<<".Easy";
            ss >> name;
            vector < double > vec;
            ret = reader.readNumberArray( name , vec );
            problem( name , ret );
            if (vec.size() == 3 ){
                double ez[3] = {vec[0] , vec[1] , vec[2] };
                s->setEasyAngles( ez );
                s->calcV1V2(); // calculates surface vectors from easy angles
            }
            else{
                cout << "error - Easy direction "<< name <<" needs 3 angles - bye!" << endl;
                exit(1);
            }

            ss.clear(); name.clear();
            ss << "FIXLC" << i << ".K1";
            ss >> name;
            ret = reader.readNumber( name , dbl_val);
            if (! s->getisFixed () ){
                problem(name , ret);
                s->setK1(dbl_val);
            }
            ss.clear(); name.clear();
            ss << "FIXLC" << i << ".K2";
            ss >> name;
            ret = reader.readNumber( name , dbl_val);
            if (! s->getisFixed () ){
                problem(name , ret);
                s->setK2(dbl_val);
            }
            alignment->addSurface(s);
        }
        // end if FXLCi exists
    }// end for every FIXLC#

}//end void readAlignment

void readElectrodes(Electrodes* electrodes,
                    Reader& reader,
                    EventList& evli)
{

    std::string name;
    std::stringstream ss;
    // COUNT NUMBER OF ELECTRODES
    int i = 1;
    for ( ; i < 100 ; i++){
        std::vector <double> dummy;
        name.clear();
        ss.clear();
        ss << "E"<< i <<".Pot";
        ss >> name;
        int ret = reader.readNumberArray(name,dummy );

        if (ret!=READER_SUCCESS){
            break;
        }
    }
    size_t numElectrodes = i - 1 ;
    electrodes->setnElectrodes( numElectrodes );

    for (size_t i = 1 ; i < numElectrodes + 1 ; i++){
        std::vector<double> times;
        std::vector<double> pots;

        name.clear();
        ss.clear();
        ss << "E"<<i<<".Time";
        ss >> name;
        int ret = reader.readNumberArray(name,times);
        problem(name, ret);

        name.clear();
        ss.clear();
        ss << "E"<<i<<".Pot";
        ss >> name;
        ret = reader.readNumberArray(name, pots);

        if (times.size() != pots.size() ){
            cout <<" error reading E" << i <<" Time/Pot lengths do not match - bye!"<<endl;
            exit(1);
        }


        // TODO: DISTINGUISH BETWEEN TIME/ITERATION SWITCHING
        for (size_t j = 0 ; j < times.size() ; j++){
            Event* swEvent = new Event(EVENT_SWITCHING, times[j]);
            SwitchingInstance* si = new SwitchingInstance(times[j],  // WHEN
                                                          pots[j],   // NEW POTENTIAL VALUE
                                                          i-1);      // WHICH ELECTRODE
            swEvent->setEventDataPtr( static_cast<void*>(si) );
            evli.insertTimeEvent( swEvent );
        }
    }

    // READ DIELECTRIC PERMITTIVITIES

    name = "eps_dielectric";
    std::vector<double> eps_temp;
    int ret = reader.readNumberArray(name, eps_temp );
    problem_format(name, ret );
    if ( ret == READER_SUCCESS ){
        electrodes->eps_dielectric = eps_temp;
    }


    // READ UNIFORM ELECTRIC FIELD
    std::vector<double> Efield;
    name = "EField";
    ret = reader.readNumberArray(name, Efield);
    problem_format( name , ret);
    if ( ( ret == READER_SUCCESS ) &&
         ( Efield.size() == 3 ) ){
        electrodes->EField[0] = Efield[0];
        electrodes->EField[1] = Efield[1];
        electrodes->EField[2] = Efield[2];

        // TODO ADD SWITCHING INSTANCE WITH SPECIAL ELECTRODE NUMBER
        // INDICATING THAT A UNIFORM ELECTRC FIELD WILL BE CONSIDERED
        Event* swEvent = new Event( EVENT_SWITCHING , 0.0 );
        SwitchingInstance* si = new SwitchingInstance( 0 ,
                                                       0 ,
                                                       SwitchingInstance::UNIFORM_E_FIELD
                                                       );
        swEvent->setEventDataPtr( static_cast <void*> (si) );
        evli.insertTimeEvent( swEvent );
    }

    electrodes->setImplicitVariables();
}
// end readElectrodes

void readMeshrefinement( MeshRefinement* meshrefinement, Reader& reader){
    int num_refreg = 100;
    stringstream ss;
    string name = "";
    std::vector < double > vec;

    std::string str_val;
    int ret = 0;

    // READS REFREG STRUCTURES
    for ( int i = 0 ; i < num_refreg ; i++){
        ss.clear(); name.clear();
        ss << "REFREG" << i << ".Type";
        ss >> name;

        ret = reader.readString(name , str_val);

        if (ret == READER_SUCCESS ){ // if this REFREG was found
            cout << " reading REFREG" << i << endl;
            RefReg rr;
            rr.setType(str_val);

            // READ REGION COORDINATES
            // READ X
            ss.clear(); name.clear();
            ss << "REFREG" << i <<".X";
            ss >> name;
            ret = reader.readNumberArray(name, vec);
            problem(name, ret);
            rr.addX( vec );
            // READ Y
            ss.clear(); name.clear();
            ss << "REFREG" << i <<".Y";
            ss >> name;
            ret = reader.readNumberArray(name, vec);
            problem(name, ret);
            rr.addY( vec );
            // READ Z
            ss.clear(); name.clear();
            ss << "REFREG" << i <<".Z";
            ss >> name;
            ret = reader.readNumberArray(name, vec);
            problem(name, ret);
            rr.addZ( vec );
            // READ DISTANCE
            ss.clear(); name.clear();
            ss << "REFREG" << i <<".Distance";
            ss >> name;
            ret = reader.readNumberArray(name, vec);
            cout << name << " " << ret << "length vec = " << vec.size() << endl;
            problem(name, ret);
            rr.addDistance( vec );
            //rr.PrintRefinementRegion();
            meshrefinement->addRefinementRegion( rr );
        }
    }
}//end readMeshrefinement
void readAutorefinement( MeshRefinement* meshrefinement, Reader& reader){
    int num_autoref = 100;
    stringstream ss;
    string name = "";
    std::vector <double> vec_val;
    std::string str_val;
    double dbl_val;
    int ret = 0;
    for (int i = 0 ; i < num_autoref; i++){ // fro autorefs, i
        ss.clear();
        name.clear();

        ss <<"AUTOREF"<< i << ".Type";
        ss >> name;

        ret = reader.readString(name, str_val);

        if (ret == READER_SUCCESS){ // IF THIS AUTOREFINEMENT IS FOUND
            cout << "reading AUTOREF"<< i << endl;
            AutoRef ar;
            if(!ar.setType(str_val)){
                cout << "error - invalid AUTOREF"<<i << " Type: " << str_val <<" - bye!" << endl;
                exit(1);
            }

            ss.clear(); name.clear();
            ss << "AUTOREF" << i <<".RefIter";
            ss >> name;
            ret = reader.readNumber(name, dbl_val );
            problem(name, ret);
            ar.setRefIter((int) dbl_val );

            // READ MaxValue
            ss.clear(); name.clear();
            ss << "AUTOREF" << i <<".MaxValue";
            ss >> name;
            ret = reader.readNumberArray( name, vec_val );
            problemo(name, ret );
            ar.addMaxValue( vec_val );

            // READ MinSize
            ss.clear(); name.clear();
            ss << "AUTOREF"<<i<<".MinSize";
            ss >> name;
            ret = reader.readNumber( name , dbl_val);
            problemo ( name, ret );
            ar.setMinSize(dbl_val);
            printf("Autoref.MinSize = %e\n", dbl_val);
            meshrefinement->addAutorefinement( ar );

        }// end if autoref i found
    }// end for possible autorefs, i

}

void readEndrefinement( MeshRefinement* meshrefinement, Reader& reader){
    int num_endref = 10;
    stringstream ss;
    string name = "";
    std::vector <double> vec_val;
    std::string str_val;
    double dbl_val;
    int ret = 0;
    for (int i = 0 ; i < num_endref; i++){ // fro autorefs, i
        ss.clear();
        name.clear();

        ss <<"ENDREF"<< i << ".Type";
        ss >> name;
        //cout << name << endl;
        ret = reader.readString(name , str_val);
        //printf(" ret = %i\n", ret );
        if (ret == READER_SUCCESS){ // IF THIS AUTOREFINEMENT IS FOUND
            cout << "reading ENDREF"<< i << endl;
            EndRef er;
            if(!er.setType(str_val)){
                cout << "error - invalid ENDREF"<<i << " Type: " << str_val <<" - bye!" << endl;
                exit(1);
            }

            ss.clear(); name.clear();
            ss << "ENDREF" << i <<".RefIter";
            ss >> name;
            ret = reader.readNumber(name, dbl_val );
            problem(name, ret);
            er.setRefIter((int) dbl_val );

            // READ MaxValue
            ss.clear(); name.clear();
            ss << "ENDREF" << i <<".MaxValue";
            ss >> name;
            ret = reader.readNumberArray( name, vec_val );
            problemo(name, ret );
            er.addMaxValue( vec_val );

            // READ MinSize
            ss.clear(); name.clear();
            ss << "ENDREF"<<i<<".MinSize";
            ss >> name;
            ret = reader.readNumber( name , dbl_val);
            problemo ( name, ret );
            er.setMinSize(dbl_val);
            printf("Endref.MinSize = %e\n", dbl_val);
            meshrefinement->addEndrefinement( er );

        }// end if autoref i found
    }// end for possible autorefs, i
}

void readRefinement(Reader& reader,
                    EventList& evli)
{
    // READS SETINGS FOR MESH REFINEMENT AND ADDS REFINEMENT EVENT TO EVENT QUEUE

    // FIRST READ "GLOBAL" REPEATED REFINMENT SETTINGS
    string key = "RepRefIter";
    size_t  repRefIter(0);
    int ret = reader.readNumber( key, repRefIter );
    problem_format( key, ret );
    evli.setRepRefIter( repRefIter );

    key = "RepRefTime";
    double repRefTime(0);
    ret = reader.readNumber( key , repRefTime );
    problem_format(key, ret );
    evli.setRepRefTime( repRefTime );

    // LOOP OVER POSSIBLE REFINEMENT SETTINGS VALUES
    unsigned int max_num_ref = 100;
    for (unsigned int i = 1 ; i <= max_num_ref ; i++){ // FOR REFINEMENTS
        string type = "";    // SETTING RETURN STRING VALUE
        key = setStructureKey("REFINEMENT", i , "Type");
        int ret = reader.readString(key, type);


        if ( ret == READER_SUCCESS ){ // IF REFINEMENTi IS DEFINED IN SETTINGS FILE
            // OPTIONAL, EXPLICIT ITERATIONS WHEN REFINEMENT OCCURS
            key = setStructureKey("REFINEMENT", i , "Iterations");
            vector <long int> iterations;
            ret = reader.readNumberArray( key, iterations );
            problem_format( key , ret );

            // OPTIONAL, EXPLICIT TIMES WHEN REFINEMENT OCCURS
            key = setStructureKey("REFINEMENT",i,"Times");
            vector <double> times;
            reader.readNumberArray( key, times );
            problem_format( key , ret );

            // MAKE SURE THAT ONLY Iterations OR Times IS DEFINED, NOT BOTH
            if ( (iterations.empty() ) &&
                 (!times.empty() ) ){
                char msg[200];
                sprintf( msg, "error with REFINEMENT%i, can't define both Iterations AND Times for same setting - bye!", i);
                problem( msg );
            }

            // START READING ALL POSSIBLE RefInfo DATA FIELDS
            key = setStructureKey("REFINEMENT", i , "Values");
            vector<double> values;
            ret = reader.readNumberArray(key, values);
            problem_format(key, ret);

            vector<double> X , Y, Z;
            key = setStructureKey("REFINEMENT", i , "X");
            ret = reader.readNumberArray( key, X );
            problem_format( key, ret );
            key = setStructureKey("REFINEMENT", i , "Y");
            ret = reader.readNumberArray( key, Y );
            problem_format( key, ret );
            key = setStructureKey("REFINEMENT", i , "Z");
            ret = reader.readNumberArray( key, Z );
            problem_format( key, ret );

            // DATA HAS BEEN COLLECTED. DEPENDING ON WHETHER EXPLICIT REFINEMENT
            // EVENTS HAVE BEEN DEFINED ADD EVENT(s) TO EVENT LIST

            // IF EXPLICIT ITERATIOSN ARE DEFINED, BREAK THEM TO SEPARATE EVENT
            if (iterations.empty() ) {
                for (size_t j = 0 ; j < iterations.size() ; j++){
                    unsigned int itr = (unsigned int) iterations[j];
                    RefInfo* refinfo = new RefInfo(type);
                    refinfo->setIteration( itr );
                    refinfo->setValues( values );
                    refinfo->setCoords(X, Y, Z);
                    RefInfo::validate( *refinfo );
                    Event* e = new Event(EVENT_REFINEMENT, itr);
                    e->setEventDataPtr( static_cast<void*> (refinfo) );
                    evli.insertIterEvent( e );
                }

            }
            // IF EXPLICIT TIMES ARE DEFINED, BREAK INTO SEPARATE EVENTS
            else if ( !times.empty() ){
                for (size_t j = 0 ; j < times.size() ; j++){
                    double tme = times[j];
                    RefInfo* refinfo = new RefInfo(type);
                    refinfo->setTime( tme );
                    refinfo->setValues( values );
                    refinfo->setCoords(X, Y, Z);
                    RefInfo::validate( *refinfo );
                    Event* e = new Event( EVENT_REFINEMENT, tme );
                    e->setEventDataPtr( static_cast<void*> (refinfo) );
                    evli.insertTimeEvent( e );
                }
            }
            // REPEATING REFINEMENT
            else{
                RefInfo* refinfo = new RefInfo(type);
                refinfo->setValues( values );
                refinfo->setCoords(X, Y, Z);
                RefInfo::validate( *refinfo );
                Event* e = new Event(EVENT_REFINEMENT, refinfo->getRefIter() ); // REPETITION BY ITERATION ONLY CURRENTLY SUPPORTED
                e->setEventDataPtr( static_cast<void*> (refinfo) );
                evli.addRepRefInfo( e );
            }

        }// END IF SUCCESS
    }// END FOR LOOP OVER REFINEMENTi
}

void ReadSettings(
        string settings_filename,
        Simu* simu,
        LC& lc,
        Boxes* boxes,
        Alignment* alignment,
        Electrodes* electrodes,
        MeshRefinement* meshrefinement,
        EventList& eventlist)
{

    Reader reader;
    reader.isIgnoreCase = true;
    using namespace std;
    if ( reader.openFile(settings_filename) ){
        cout << "reading: "<< settings_filename << endl;
        // READ SIMU
        readSimu(simu , reader, eventlist);

        // READ LC
        readLC(lc , reader);
        lc.convert_params_n2Q();

        // READ BOXES
        readBoxes(boxes , reader);

        // READ ALIGNMENT SURFACES
        readAlignment(alignment, reader);

        // READ ELECTRODES
        readElectrodes(electrodes , reader, eventlist);

        // READ REFINEMENT
        readRefinement( reader, eventlist );

        // exit(0);

        reader.closeFile();
    }
    else{
        cout << "error - could not open " << settings_filename <<" - bye!" << endl;
        exit(1);
    }
}
// end ReadSettings



void ReadSolverSettings(const char* filename, Settings* settings)
{
    // READS SOLVER SETTINGS FROM FILE. THIS WAS ORIGIANLLY ASSUMED TO BE
    // CALLED "solver.qfg", BUT NOW ALSO READS "solver.txt", IF "solver.qfg"
    // IS NOT FOUND
    std::string fn = filename;
    if ( !FilesysFun::fileExists(fn) )
    {
        fn = "solver.txt";
    }


    Reader reader;
    if ( reader.openFile(fn) )
    {

        int i_val;
        double dbl;
        int ret;
        string name;

        name = "nThreads";
        ret = reader.readNumber(name, i_val);
        problem(name, ret);
        settings->setnThreads(i_val);

        name = "Q_Solver";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_Solver( i_val);

        name = "V_Solver";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setV_Solver( i_val);

        name = "Q_Newton_Panic_Iter";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_Newton_Panic_Iter( i_val);

        name = "Q_Newton_Panic_Coeff";
        problem(name , reader.readNumber(name , dbl ) );
        settings->setQ_Newton_Panic_Coeff(dbl);


        name = "Q_PCG_Preconditioner";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_PCG_Preconditioner( i_val);

        name = "Q_PCG_Maxiter";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_PCG_Maxiter( i_val);


        name = "Q_PCG_Toler";
        problem(name , reader.readNumber(name , dbl ) );
        settings->setQ_PCG_Toler(dbl);

        name = "Q_GMRES_Preconditioner";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_GMRES_Preconditioner(i_val);


        name = "Q_GMRES_Maxiter";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_GMRES_Maxiter( i_val);


        name = "Q_GMRES_Restart";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setQ_GMRES_Restart(i_val);

        name = "Q_GMRES_Toler";
        problem(name , reader.readNumber(name , dbl) );
        settings->setQ_GMRES_Toler(dbl);

        name = "V_PCG_Preconditioner";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setV_PCG_Preconditioner( i_val);

        name = "V_PCG_Maxiter";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setV_PCG_Maxiter(i_val);

        name = "V_PCG_Toler";
        problem(name , reader.readNumber(name , dbl ) );
        settings->setV_GMRES_Toler(dbl);

        name = "V_GMRES_Preconditioner";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setV_GMRES_Preconditioner( i_val);

        name = "V_GMRES_Maxiter";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setV_GMRES_Maxiter( i_val);

        name = "V_GMRES_Restart";
        problem(name , reader.readNumber(name , i_val ) );
        settings->setV_GMRES_Restart( i_val);

        name = "V_GMRES_Toler";
        problem(name , reader.readNumber(name , dbl ) );
        settings->setV_GMRES_Toler( dbl );

        //settings->PrintSettings();

        reader.closeFile();

    }
    else
    {
        cout << "Did not read solver settings. Using defaults" << endl;
        //cout<< "Could not read solver settings file:" << fn << " - bye!"<< endl;
        //exit(1);
    }

} // end readSolverSettings

