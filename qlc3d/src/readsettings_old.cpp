#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h>
#include <qlc3d.h>
#include <string>
#include <vector>
#include <reader.h>
#include <iostream>
using namespace std;


void problem(const char *error)
{

	printf("WARNING - problem reading setting for : ");
	printf("%s\n",error);
	//printf("\n");
	//exit(1);
}
void problem(std::string& name, int ret){
    Reader temp;
    if (ret!= READER_SUCCESS){
    cout << "Problem reading " << name << " , computer says: " << temp.getErrorString(ret) << endl;
    }
}
void problemo(std::string& name, int ret){
    Reader temp;
    if ( (ret!= READER_SUCCESS) && (ret!=READER_NOT_FOUND)){
    cout << "Problem reading " << name << " , computer says: " << temp.getErrorString(ret) << endl;
    }
}
// end problemo (optional parameter problem)


void readLC(LC* lc, config_t* cfg)
{
	config_setting_t *setting;
	config_set_auto_convert(cfg,1);

	// ELASTIC COEFFS
			setting = config_lookup(cfg,"LC.K11");
				if (setting) lc->K11 = config_setting_get_float(setting);
				else problem("LC.K11");
			setting = config_lookup(cfg,"LC.K22");
				if (setting) lc->K22 = config_setting_get_float(setting);
				else problem("LC.K22");
			setting = config_lookup(cfg,"LC.K33");
				if (setting) lc->K33 = config_setting_get_float(setting);
				else problem("LC.K33");
			setting = config_lookup(cfg,"LC.p0");
				if (setting) lc->p0 = config_setting_get_float(setting);
				else problem("LC.p0");

	// THERMOTROPIC COEFFS
			setting = config_lookup(cfg,"LC.A");
				if (setting) lc->A = config_setting_get_float(setting);
				else problem("LC.A");
			setting = config_lookup(cfg,"LC.B");
				if (setting) lc->B = config_setting_get_float(setting);
				else problem("LC.B");
			setting = config_lookup(cfg,"LC.C");
				if (setting) lc->C = config_setting_get_float(setting);
				else problem("LC.C");
	//ELECTRIC COEFFICIENTS
			setting = config_lookup(cfg,"LC.eps_par");
				if (setting) lc->eps_par = config_setting_get_float(setting);
				else problem("LC.eps_par");
			setting = config_lookup(cfg,"LC.eps_per");
				if (setting) lc->eps_per = config_setting_get_float(setting);
				else problem("LC.eps_per");
			setting = config_lookup(cfg,"LC.e11");
				if (setting) lc->e11 = config_setting_get_float(setting);
				else problem("LC.e11");
			setting = config_lookup(cfg,"LC.e33");
				if (setting) lc->e33 = config_setting_get_float(setting);
				else problem("LC.e33");
	//VISCOUS COEFFICIENTS
			setting = config_lookup(cfg,"LC.gamma1");
				if (setting) lc->gamma1 = config_setting_get_float(setting);
				else problem("LC.gamma1");
			setting = config_lookup(cfg,"LC.gamma2");
				if (setting) lc->gamma2 = config_setting_get_float(setting);
				else problem("LC.gamma2");
			setting = config_lookup(cfg,"LC.alpha1");
				if (setting) lc->alpha1 = config_setting_get_float(setting);
				else problem("LC.alpha1");
			setting = config_lookup(cfg,"LC.alpha4");
				if (setting) lc->alpha4 = config_setting_get_float(setting);
				else problem("LC.alpha4");
			setting = config_lookup(cfg,"LC.alpha5");
				if (setting) lc->alpha5 = config_setting_get_float(setting);
				else problem("LC.alpha5");
			setting = config_lookup(cfg,"LC.alpha6");
				if (setting) lc->alpha6 = config_setting_get_float(setting);
				else problem("LC.alpha6");

}//end void readLC


void readSimu(Simu* simu, Reader& reader)
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
    problemo(name, ret);

    name = "SaveDir";
    ret = reader.readString(name , str_var);
    if ( ret == READER_SUCCESS)
        simu->setSaveDir( str_var);
    problemo(name , ret);
//========================
// SCALAR VALUES
//========================
    name = "EndValue";
    ret = reader.readNumber(name, dbl_var);
    if (ret == READER_SUCCESS)
        simu->setEndValue(dbl_var);
    problem(name, ret);

    name = "dt";
    ret = reader.readNumber(name , dbl_var);
    if (ret == READER_SUCCESS)
        simu->setdt(dbl_var);
    problemo(name, ret);

    name = "Maxdt";
    ret = reader.readNumber(name , dbl_var);
    if(ret==READER_SUCCESS)
        simu->setMaxdt(dbl_var);
    problemo(name, ret);

    name = "MaxError";
    ret = reader.readNumber(name, dbl_var);
    if(ret==READER_SUCCESS)
        simu->setMaxError(dbl_var);
    problemo(name, ret);

    name = "OutputEnergy";
    ret  = reader.readNumber(name , int_var);
    if ( ret == READER_SUCCESS)
        simu->setOutputEnergy(int_var);
    problemo(name, ret);

//------------------------------------------
//
// READ VECTOR VALUES
//
//------------------------------------------
    name = "StretchVector";
    vector <double> vec;
    ret = reader.readNumberArray(name , vec );
    if (ret == READER_SUCCESS){
        simu->setStretchVectorX(vec[0]);
        simu->setStretchVectorY(vec[1]);
        simu->setStretchVectorZ(vec[3]);
    }
    problemo(name, ret);

    simu->PrintSimu();

	/*
	const char *str;
	double 	dval;
	int 	ival;

	config_lookup_string(cfg,"SIMU.MeshName",&str);
	simu->MeshName = str;

	config_lookup_float(cfg,"SIMU.dt",&dval);
	simu->dt = dval;

	config_lookup_float(cfg,"SIMU.Maxdt",&dval);
	simu->setMaxdt(dval);

	config_lookup_float(cfg,"SIMU.MaxError",&dval);
	simu->setMaxError(dval);

	config_lookup_string(cfg,"SIMU.EndCriterion",&str);
	simu->setEndCriterion(str);
	if ( str[0] == 'I' ) // if iterations
		{
		config_lookup_int(cfg,"SIMU.EndValue",&ival);
		simu->setEndValue(ival);
		}
	else
	if (( str[0] == 'T' ) || (str[0] == 'C')) // if Time or Change
		{
		config_lookup_float(cfg,"SIMU.EndValue",&dval);
		simu->setEndValue(dval);
		}
	if (config_lookup_string(cfg,"SIMU.LoadQ", &str)){
		string lqs = str;
		simu->setLoadQ(lqs);
	}

// READ OPTIONAL INTEGER VALUED SETTINGS - VALID DEFAULTS HAVE BEEN SET IN CONSTRUCTOR
	if (config_lookup_int(cfg,"SIMU.OutputEnergy", &ival)){
			simu->setOutputEnergy(ival);
	}
	if (config_lookup_int(cfg,"SIMU.OutputFormat", &ival) ){
		simu->setOutputFormat(ival);
	}


// READ STRETCHVECTOR
	string setting_name = "SIMU.StretchVector";
	config_setting_t *setting;
	setting = config_lookup(cfg,setting_name.c_str());
	if (setting){
		int c1 = config_setting_length(setting);	// make sure 3 vector-components are defined
		if (c1!=3){
			printf("error - StretchVector must contain 3 components - bye \n");
			exit(1);
		}

		simu->setStretchVectorX( config_setting_get_float_elem(setting,0) );
		simu->setStretchVectorY( config_setting_get_float_elem(setting,1) );
		simu->setStretchVectorZ( config_setting_get_float_elem(setting,2) );
	}

	// READ OPTIONAL "SaveDir" setting
	setting_name = "SIMU.SaveDir";
	setting = config_lookup( cfg,setting_name.c_str() );
	string s_var;
	if (setting) // if optional setting is defined
	{
		s_var = config_setting_get_string( setting);
		s_var = s_var + "/";
		simu->setSaveDir( s_var );
		//printf("SaveDir is defined : %s \n", simu->getSaveDir().c_str());
	}

	// READ OPTIONAL "LoadDir" setting
	setting_name = "SIMU.LoadDir";
	setting = config_lookup( cfg, setting_name.c_str() );

	if (setting) // if optional setting is defined in file
	{
		s_var = config_setting_get_string( setting );
		s_var = s_var + "/";
		simu->setLoadDir( s_var );
		printf("LoadDir is defined : %s \n", simu->getLoadDir().c_str() );
	}


*/


}//end void readSimu

void readBoxes(Boxes* boxes, config_t* cfg)
{
	config_setting_t *setting;

	//limit number of boxes to 0-#0000FF99
	int num_boxes = 99;


	std::string box_name;
	std::string box_setting_name;
	char number[3];


	for (int i = 0; i<= num_boxes ; i++) // loop over each box
		{
			//itoa(i,number,10); // convert number to text
			sprintf(number,"%i",i);
			box_name = "LC_ORIENTATION.BOX";// + number;
			box_name.append(number);

			setting = config_lookup(cfg,box_name.c_str());
			if (setting) // if setting found
			{

				//printf("Reading ");
				//printf(box_name.c_str());
				//printf(":\n");

				boxes->addBox();

				box_setting_name = box_name + ".Type";
				setting= config_lookup(cfg,box_setting_name.c_str());
					if (setting) boxes->box[boxes->n_Boxes-1]->Type = config_setting_get_string(setting);
					else problem("BOX1.Type");
				//*
				box_setting_name =box_name+".X";
				setting = config_lookup(cfg,box_setting_name.c_str());
					if (setting)
					{
					boxes->box[boxes->n_Boxes-1]->X[0] = config_setting_get_float_elem(setting,0);
					boxes->box[boxes->n_Boxes-1]->X[1] = config_setting_get_float_elem(setting,1);
					}

				box_setting_name =box_name+".Y";
				setting = config_lookup(cfg,box_setting_name.c_str());
					if (setting)
					{
					boxes->box[boxes->n_Boxes-1]->Y[0] = config_setting_get_float_elem(setting,0);
					boxes->box[boxes->n_Boxes-1]->Y[1] = config_setting_get_float_elem(setting,1);
					}

				box_setting_name =box_name+".Z";
				setting = config_lookup(cfg,box_setting_name.c_str());
					if (setting)
					{
					boxes->box[boxes->n_Boxes-1]->Z[0] = config_setting_get_float_elem(setting,0);
					boxes->box[boxes->n_Boxes-1]->Z[1] = config_setting_get_float_elem(setting,1);
					}

				box_setting_name =box_name+".Tilt";
				setting = config_lookup(cfg,box_setting_name.c_str());
					if (setting)
					{
					boxes->box[boxes->n_Boxes-1]->Tilt[0] = config_setting_get_float_elem(setting,0);
					boxes->box[boxes->n_Boxes-1]->Tilt[1] = config_setting_get_float_elem(setting,1);
					}
				box_setting_name =box_name+".Twist";
				setting = config_lookup(cfg,box_setting_name.c_str());
					if (setting)
					{
					boxes->box[boxes->n_Boxes-1]->Twist[0] = config_setting_get_float_elem(setting,0);
					boxes->box[boxes->n_Boxes-1]->Twist[1] = config_setting_get_float_elem(setting,1);
					}
				//
			}
				else
				{
					//do nothing
					//printf("not found\n");
				}


		}// end for i



}// end void readBoxes

void readAlignment(Alignment* alignment, config_t* cfg)
{
	config_setting_t *setting;
	int num_surfaces = 99; // limit surfaces to 0-99
	std::string surface_name;
	std::string surface_setting_name;
	char number[3];


	for ( int i = 0; i < num_surfaces; i++)
	{
			//itoa(i,number,10); // convert number to text
			sprintf(number,"%i",i);
			surface_name = "ALIGNMENT.FIXLC";// + number;
			surface_name.append(number);

			setting = config_lookup(cfg,surface_name.c_str());
			if (setting) // if setting found
			{
				const char *str;
				double 		dval;
				alignment->addSurface();

				//READ ANCHORING TYPE
				surface_setting_name = surface_name;
				surface_setting_name.append(".Anchoring");
				config_lookup_string(cfg,surface_setting_name.c_str(),&str);
				alignment->surface[alignment->getnSurfaces()-1]->setAnchoringType( str);

				//READ ANCHORING STRENGTH
				surface_setting_name = surface_name;
				surface_setting_name.append(".Strength");
					config_lookup_float(cfg,surface_setting_name.c_str(),&dval);
				alignment->surface[alignment->getnSurfaces()-1]->setStrength( dval );

				//READ K1
				surface_setting_name = surface_name;
				surface_setting_name.append(".K1");
				config_lookup_float(cfg,surface_setting_name.c_str(),&dval);
				alignment->surface[alignment->getnSurfaces()-1]->setK1(dval);

				//READ K2
				surface_setting_name = surface_name;
				surface_setting_name.append(".K2");
				config_lookup_float(cfg,surface_setting_name.c_str(),&dval);
				alignment->surface[alignment->getnSurfaces()-1]->setK2(dval);

				//READ EASY DIRECTION ANGLES
				surface_setting_name = surface_name;
				surface_setting_name.append(".Easy");
				setting = config_lookup(cfg,surface_setting_name.c_str());
				if (setting)
				{
					double temp[3];
					temp[0] = config_setting_get_float_elem(setting,0);
					temp[1] = config_setting_get_float_elem(setting,1);
					temp[2] = config_setting_get_float_elem(setting,2);
					alignment->surface[alignment->getnSurfaces()-1]->setEasyAngles(temp);
				}
			}



	}

}//end void readAlignment

void readElectrodes(Electrodes* electrodes, config_t* cfg)
{
	config_setting_t *setting;
	int num_electrodes = 99; // limit electrodes to 0-99
	std::string electrode_name;	// E1, E2, E3 ...
	std::string electrode_setting_name; // either Time or Pot
	char number[3];


	for ( int i = 0; i < num_electrodes; i++)	// read electrodes
	{

			sprintf(number,"%i",i);
			electrode_name = "ELECTRODES.E";// + number;
			electrode_name.append(number);

			setting = config_lookup(cfg,electrode_name.c_str());

			if (setting) // if setting found
			{
					electrode_setting_name = electrode_name;
					electrodes->AddElectrode();
					int current_electrode_num = electrodes->getnElectrodes()-1;




				// read switching time array
				int c1= 0; // number of switching times
				int c2= 0; // number of switching potentials
				electrode_setting_name = electrode_name;
				electrode_setting_name.append(".Time");
				setting = config_lookup(cfg,electrode_setting_name.c_str());

				if (setting)
				{

					c1 = config_setting_length(setting);	// count number of switching times

					for (int j = 0; j <c1; j++)
					{
						electrodes->E[current_electrode_num]->Time.push_back( config_setting_get_float_elem(setting,j));
					}
				}
				else
				{
					printf("problem with setting file - could not find times for electrode %i\n",i);
					exit(1);
				}

				// read potential array
				electrode_setting_name = electrode_name;
				electrode_setting_name.append(".Pot");
				setting = config_lookup(cfg,electrode_setting_name.c_str());
				if (setting)
				{
					c2 = config_setting_length(setting);	// count number of switching times
					if (c2!=c1)
					{
						printf("number of switching times must equal number of potentials, bye!");
						exit(1);
					}
					for (int j = 0; j <c2; j++)
					{
						electrodes->E[electrodes->nElectrodes-1]->Potential.push_back(config_setting_get_float_elem(setting,j));
					}

				}
				else
				{
					printf("problem with setting file - could not find potentials for electrode %i\n",i);
					exit(1);
				}

			}
	}
	// read dielectric permittivity values
	std::string eps =  "ELECTRODES.eps_dielectric";
	setting = config_lookup(cfg,eps.c_str());

	if (setting)// if permittivities defined
	{
		int c = config_setting_length(setting); // count number of dielectric materials
		for (int i = 0 ; i < c ; i++ )
		{
			electrodes->eps_dielectric.push_back( config_setting_get_float_elem(setting, i) );	// stores dielectric values
		}
	}


}

void readMeshRefinement(MeshRefinement* meshrefinement, config_t* cfg)
{
	config_setting_t *setting;
	int num_regions = 99; // limit refinment regions to 0-99
	std::string region_name;
	std::string region_setting_name;
	char number[3];


	for ( int i = 0; i < num_regions; i++)
	{

			//itoa(i,number,10); // convert number to text
			sprintf(number,"%i",i);
			region_name = "MESHREFINEMENT.REGION";// + number;
			region_name.append(number);



			setting = config_lookup(cfg,region_name.c_str());
			if (setting) // if setting found
			{
				//printf("reading %s\n", region_name.c_str());
				const char *str;
				//double 		dval;
				meshrefinement->addRefinementRegion();

				vector <RefReg>::iterator ritr;
				ritr = meshrefinement->RefinementRegion.end()-1; // pointer to last region

			//READ REFINEMENT TYPE
				region_setting_name = region_name;
				region_setting_name.append(".Type");
				config_lookup_string(cfg,region_setting_name.c_str(),&str);

				string s = str;
				ritr->setType( s );

			//READ PARAMS
				region_setting_name = region_name;
				region_setting_name.append(".Params");
				setting = config_lookup(cfg,region_setting_name.c_str() );
				if (setting)
				{
					int num_params = config_setting_length(setting) ; // number of params

					for (int n = 0 ; n < num_params ; n ++)
						ritr->addParams(config_setting_get_float_elem(setting,n));
				}

			// READ DISTANCES
					int num_distances;
					region_setting_name = region_name;
					region_setting_name.append(".Distance");
					setting = config_lookup(cfg,region_setting_name.c_str() );
					if(setting)
					{
						num_distances = config_setting_length(setting);
						for (int n = 0 ; n < num_distances ; n++ )
							ritr->addDistance(config_setting_get_float_elem(setting,n));
					}
					else
					{
						printf("error -%s\nDistance is not defined - bye!\n",region_name.c_str());
						exit(1);
					}
			// READ COORDINATES
				int num_x = 0;
				int num_y = 0;
				int num_z = 0;
			// X
				region_setting_name = region_name;
				region_setting_name.append(".X");
				setting = config_lookup(cfg , region_setting_name.c_str() );
				if (setting)
				{
					num_x = config_setting_length(setting); // get number x - coordinates
					for(int n = 0 ; n < num_x ; n ++ )
						ritr->addX(config_setting_get_float_elem(setting,n));
				}
			// Y
				region_setting_name = region_name;
				region_setting_name.append(".Y");
				setting = config_lookup(cfg , region_setting_name.c_str() );
				if (setting)
				{
					num_y = config_setting_length(setting); // get number y - coordinates
					for(int n = 0 ; n < num_y ; n ++ )
						ritr->addY(config_setting_get_float_elem(setting,n));
				}
			// Z
				region_setting_name = region_name;
				region_setting_name.append(".Z");
				setting = config_lookup(cfg , region_setting_name.c_str() );
				if (setting)
				{
					num_z = config_setting_length(setting); // get number z - coordinates
					for(int n = 0 ; n < num_z ; n ++ )
						ritr->addZ(config_setting_get_float_elem(setting,n));
				}
			if ( (num_x!=num_y) && (num_x != num_z) )
			{
				printf("error - %s\nnumber of x,y,z coordinates must be same - bye!\n",region_name.c_str());
				exit(1);
			}


			}



	}

}//end void readMeshRefinement

void ReadSettings(string settings_filename, Simu* simu, LC* lc, Boxes* boxes, Alignment* alignment, Electrodes* electrodes, MeshRefinement* meshrefinement)
{


//struct config_t cfg;
//struct config_setting_t *setting;
//FILE *fid;

//fid= fopen(settings_filename.c_str() , "rt" );
//if (fid==NULL)
//{
//	printf("could not open settings file: %s - bye!\n" , settings_filename.c_str());
//	exit(1);
//}
//config_init(&cfg);

    Reader reader;
    using namespace std;
    if ( reader.openFile(settings_filename) ){
        cout << "reading: "<< settings_filename << endl;
    // READ SIMU
        readSimu(simu , reader);
        reader.closeFile();
    }
    else{
        cout << "error 0 could not open " << settings_filename <<" -bye!" << endl;
        exit(1);
    }




/*
if (config_read(&cfg,fid))
	{
	printf("config file %s seems ok %c\n", settings_filename.c_str() , 001  );


	//READ SIMU CONFIGURATION
		setting = config_lookup(&cfg,"SIMU");
		if (config_lookup(&cfg,"SIMU"))
		{
			readSimu(simu,&cfg);
		}
		else
		{
			printf("not reading SIMU ?\n");
		}

	//READ LC MATERIAL PARAMETERS
		if (! config_lookup(&cfg,"LC"))
			printf("not reading LC ?\n");
		else
		{
			readLC(lc, &cfg);
			lc->convert_params_n2Q();

		}// end read LC material parameters


	// READ INITIAL CONFIGURATION BOXES
	// int num_boxes=0;
	 if(!config_lookup(&cfg,"LC_ORIENTATION"))
	 {
		printf("no initial LC orientation specified\n");
	 }
	 {
		readBoxes(boxes,&cfg);

	 }// end read boxes



	//READ ALIGNEMNT CONFIGURATION
	if(config_lookup(&cfg,"ALIGNMENT"))
	{
		readAlignment(alignment, &cfg);
	}
	else
	{
		printf("not reading ALIGNMENT ?\n");
	}


	// READ ELECTRODES CONFIGURATION
	if(config_lookup(&cfg,"ELECTRODES"))
	{
		readElectrodes(electrodes, &cfg);
	}
	else
	{
		printf("not reading ELECTRODES ?\n");
	}

	// READ MESHREFINEMENT
	if (config_lookup(&cfg,"MESHREFINEMENT"))
	{
		readMeshRefinement(meshrefinement, &cfg);
	}
	else
	{
		printf("not reading MESHREFINEMENT\n");
	}


	}
	else
	{
		//const char* cc = config_error_text(&cfg);
		printf("problem reading configuration file line %i\n", config_error_line(&cfg));
		printf("%s\n",config_error_text(&cfg) );
		exit(1);
	}// end reading settings


fclose(fid);
//return 0;
*/
}


void settingserror(int cfg_ok, string qfg)
{
	if (cfg_ok == CONFIG_FALSE)
	{
		printf("error in reading settings.qfg at %s - bye!\n",qfg.c_str());
		exit(1);
	}
}

void ReadSolverSettings(const char* filename, Settings* settings)
{
	/*
	printf("Reading solver settings from: %s\n", filename);

	FILE* fid= fopen(filename,"rt");
	if (fid==NULL)
	{
		printf("Could not open the file: %s - bye!\n", filename);
		exit(1);
	}

	struct config_t cfg;
	//struct config_setting_t* S;
	config_init(&cfg);

	if ( ! config_read(&cfg , fid ) )
	{
		printf("Something seems to be wrong with the file ?! - bye\n");
		exit(1);
	}

		if (config_lookup(&cfg,"Settings"))
		{
		int qfg_ok;
		string qfg;
		int ival;
		double   fval;
        qfg ="Settings.nThreads";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival); 		//	number of threads
		settingserror(qfg_ok,qfg);
		settings->setnThreads(ival);
		//printf("number of threads = %li\n", ival );


	qfg ="Settings.Q_Solver";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);		// 	Choose Q-tensor solver
		//printf("qfg_ok = %i , Q_Solver = %li\n", qfg_ok , ival);
		settingserror(qfg_ok,qfg);
		settings->setQ_Solver(ival);

	qfg ="Settings.V_Solver";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);
		settingserror(qfg_ok,qfg);
		settings->setV_Solver(ival);

// read Q PCG settings

	qfg ="Settings.Q_Newton_Panic_Iter";					// Newton Panic Iteration
		qfg_ok = config_lookup_int(&cfg,qfg.c_str(),&ival);
		settingserror(qfg_ok,qfg);
		settings->setQ_Newton_Panic_Iter(ival);
	qfg ="Settings.Q_Newton_Panic_Coeff";
		qfg_ok = config_lookup_float(&cfg,qfg.c_str(),&fval);
		settingserror(qfg_ok,qfg);
		settings->setQ_Newton_Panic_Coeff(fval);

	qfg ="Settings.Q_PCG_Maxiter";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);	// 	Maximum Q-PCG Iterations
		settingserror(qfg_ok,qfg);
		settings->setQ_PCG_Maxiter(ival);

	qfg ="Settings.Q_PCG_Preconditioner";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);	// Choose Q-PCG preconditioner
		settingserror(qfg_ok,qfg);
		settings->setQ_PCG_Preconditioner(ival);

	qfg ="Settings.Q_PCG_Toler"	;
		qfg_ok =config_lookup_float(&cfg,qfg.c_str(),&fval);	// Q-PCG tolerance
		settingserror(qfg_ok,qfg);
		settings->setQ_PCG_Toler(fval);

// read Q GMRES settings
	qfg ="Settings.Q_GMRES_Maxiter";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(), &ival); // Maximum GMRES iterations
		settingserror(qfg_ok,qfg);
		settings->setQ_GMRES_Maxiter(ival);

	qfg ="Settings.Q_GMRES_Preconditioner";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival); // Choose Q-GMRES preconditioner
		settingserror(qfg_ok,qfg);
		settings->setQ_GMRES_Preconditioner(ival);

	qfg ="Settings.Q_GMRES_Restart"	;
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);	// Set Q-GMRES restart iterations
		settingserror(qfg_ok,qfg);
		settings->setQ_GMRES_Restart(ival);

	qfg ="Settings.Q_GMRES_Toler";
		config_lookup_float(&cfg,qfg.c_str(),&fval);	// Q-GMRES tolerance
		settingserror(qfg_ok,qfg);
		settings->setQ_GMRES_Toler(fval);

// read V PCG settings
	qfg = "Settings.V_PCG_Maxiter"	;
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);	// 	Maximum V-PCG Iterations
		settingserror(qfg_ok,qfg);
		settings->setV_PCG_Maxiter(ival);

	qfg = "Settings.V_PCG_Preconditioner";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival);	// Choose V-PCG preconditioner
		settingserror(qfg_ok,qfg);
		settings->setV_PCG_Preconditioner(ival);

	qfg = "Settings.V_PCG_Toler";
		qfg_ok =config_lookup_float(&cfg,qfg.c_str(),&fval);	// V-PCG tolerance
		settingserror(qfg_ok,qfg);
		settings->setV_PCG_Toler(fval);

// read V GMRES settings
	qfg="Settings.V_GMRES_Maxiter";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(), &ival); // Maximum V- GMRES iterations
		settingserror(qfg_ok,qfg);
		settings->setV_GMRES_Maxiter(ival);

	qfg="Settings.V_GMRES_Preconditioner";
		qfg_ok =config_lookup_int(&cfg,qfg.c_str(),&ival); // Choose V-GMRES preconditioner
		settingserror(qfg_ok,qfg);
		settings->setV_GMRES_Preconditioner(ival);

	qfg = "Settings.V_GMRES_Restart";
		qfg_ok = config_lookup_int(&cfg,qfg.c_str(),&ival);	// Set V-GMRES restart iterations
		settingserror(qfg_ok,qfg);
		settings->setV_GMRES_Restart(ival);

	qfg = "Settings.V_GMRES_Toler";
		qfg_ok = config_lookup_float(&cfg,qfg.c_str(),&fval);	// V-GMRES tolerance
		settingserror(qfg_ok,qfg);
		settings->setV_GMRES_Toler(fval);



		}
		else
		{
			printf("Could not read Settings: - bye \n");
		}

	fclose(fid);
*/
}
