#include <meshrefinement.h>
#include <float.h>
#include <sstream>

RefReg::RefReg(){
	Type = RefReg_Sphere;
	Params.clear();
	Distance.clear();
        this->material_num = 0;
        X.clear();
	Y.clear();
	Z.clear();


}
RefReg::~RefReg()
{

}

bool RefReg::setType(const std::string &type){

    // test for valid refinement region types. return true if valid, else return false
    if (type.compare("Sphere") == 0 ){
	Type = RefReg_Sphere;
	return true;
    }
    if (type.compare("Line") == 0){
	Type = RefReg_Line;
	return true;
    }
    if (type.compare("Box") == 0){
	Type = RefReg_Box;
	return true;
    }
    if (type.find("Surface") == 0){ // if type starts with "Surface"
        Type = RefReg_Surface;

        // determine material type
        // Check if FIXLC
        if (type.find("FIXLC") < std::string::npos ){
            //Yay! FIXLC REFIENEMENT SPECIFIED
            // Determine FIXLC NUMBER
            size_t pos = type.find_first_of("123456789");
            if ( pos < std::string::npos ){

                std::stringstream ss; // convert string to number
                ss << type.at(pos);
                int fixlcnum;
                ss >> fixlcnum;
                this->material_num += FIXLCN_TO_MATNUM( fixlcnum ); // get actual material number
                //printf(" REFINEMENT NEAR FIXLC%i -> matnum = %i\n", fixlcnum, material_num);
                //exit(1);
            }
            else{
                printf("error - %s is not a REFREG.Type - bye\n", type.c_str() );
                exit(1);
            }

        }
        // Check if ELECTRODE
        else
        if (type.find("E") < std::string::npos){
            printf("Refinement near electrodes not implemented yet - bye\n");
            exit(1);
        }
        else{
            printf("error - %s is not a valid REFREG.Type - bye!\n", type.c_str() );
            exit(1);
        }

        return true;
    }
    cout << "error - " << type << "is not a valid REFREG.Type - bye!"<< endl;
    exit(1);
    return false;
}

void RefReg::addDistance(const std::vector<double> D){
    Distance.insert( Distance.end(), D.begin(), D.end() );
}
void RefReg::addX(const vector <double>& x){
    X.insert(X.end() , x.begin() , x.end() );
}
void RefReg::addY(const vector <double>& y){
    Y.insert(Y.end() , y.begin() , y.end() );
}
void RefReg::addZ(const vector <double>& z){
    Z.insert(Z.end() , z.begin() , z.end() );
}

void RefReg::PrintRefinementRegion()
{
	//printf("Type : %s\n",Type.c_str());
	
	
	vector <double> :: iterator litr;
	printf("Params:");
	for (litr = Params.begin() ; litr != Params.end() ; ++litr)
		printf(" %f",*litr);


	printf("\nDistance:");
	for (litr = Distance.begin() ; litr != Distance.end() ; ++litr)
		printf(" %f",*litr);
	
	printf("\nX:");
	for (litr = X.begin() ; litr != X.end() ; ++litr)
		printf(" %f",*litr);
		
	printf("\nY:");
	for (litr = Y.begin() ; litr != Y.end() ; ++litr)
		printf(" %f",*litr);
		
	printf("\nZ:");
	for (litr = Z.begin() ; litr != Z.end() ; ++litr)
		printf(" %f",*litr);
	printf("\n");
	
}

int RefReg::getNumIterations(){
    if ( Type != RefReg_Box )
	return (int) Distance.size();
    else
	return 1;
}
double RefReg::getDistance(const int &iter){
    if ( ( Type != RefReg_Box ) && (Distance.size() > (size_t) iter) ) // Distance makes no sense for a box
        return Distance[iter];
    else
        return 0.0;
}
//====================================
//
//  AUTOREFINEMENT CLASS
//
//====================================
AutoRef::AutoRef(){
    MinSize = 0;
    RefIter = 0;
}

unsigned int AutoRef::getNumIterations(){
/*! returns length of MaxValue setting. This is the maximum number of refinement iterations that can be performed*/
    return (int) MaxValue.size();
}
bool AutoRef::setType(const std::string &type){
    if (type.compare("Change") == 0 ){ // add types if needed
        Type = Change;
        return true;
    }
    return false;
}

void AutoRef::addMaxValue(const std::vector<double> &vals){
    MaxValue.insert( MaxValue.end(), vals.begin() , vals.end() );
}
double AutoRef::getMaxValue(const unsigned int &iter){
    if (iter<MaxValue.size() )
        return MaxValue[iter];
    else
        return DBL_MAX;

}

void AutoRef::printAutoref() {
    printf("Autoref.Type = ");
    if (this->Type == Change){
           printf("Change\n");
    }
    printf("Autoref.RefIter = %i\n", this->RefIter);

    printf("Autoref.MaxValue = ");
    for (unsigned int i = 0 ; i < this->getNumIterations() ; i++){
        printf(" %f ", this->getMaxValue( i ) );
    }

    printf("Autoref.MinSize = %f\n", this->getMinSize() );

}
//===================================
//
//  END REFINEMENT CLASS
//
//===================================
EndRef::EndRef(): AutoRef()
{
    EndRefIteration = 0;
}
EndRef::~EndRef()
{
    //~AutoRef();
}
//===============================
//
//  MESHREFINEMENT CLASS
//
//===============================


void MeshRefinement::addRefinementRegion(const RefReg reg){
	RefinementRegion.push_back(reg);
}

void MeshRefinement::addAutorefinement( AutoRef reg){
    AutoRefinement.push_back( reg );
}

void MeshRefinement::addEndrefinement( EndRef reg){
    reg.printAutoref();
    EndRefinement.push_back( reg );

}

void MeshRefinement::PrintRefinementRegions(){
	vector <RefReg>::iterator ritr;
	for (ritr = RefinementRegion.begin() ; ritr!=RefinementRegion.end() ; ++ritr)
		ritr->PrintRefinementRegion();
}

MeshRefinement::MeshRefinement(){
	RefinementRegion.clear();
    refiter = 0;
    needs_new_mesh = true;
}

MeshRefinement::~MeshRefinement(){
}

int MeshRefinement::getMaxNumRefIterations(){
    vector <RefReg>::iterator ritr;
    int max = 0;
    for (ritr = RefinementRegion.begin() ; ritr !=RefinementRegion.end() ; ++ritr)
	max = ritr->getNumIterations() > max ? ritr->getNumIterations() : max;

    return max;
}

int MeshRefinement::getMaxNumAutoRefIterations(){
/*! returns maximum number of Auto-refinement iterations */

    vector <AutoRef>::iterator itr;
    int max = 0;

    for ( itr = AutoRefinement.begin() ; itr!= AutoRefinement.end() ; ++itr){
       // itr->printAutoref();
        max = (int) itr->getNumIterations() > max ? itr->getNumIterations() : max;
       // printf(" Meshrefinement::getmaxnumautoref... ax  = %i \n", max);
    }
    return max;
}

int MeshRefinement::getMaxNumEndRefIterations(){
/*! returns maximum number of End-refinement iterations */
    vector <EndRef> ::iterator itr;
    int max = 0;
    for (itr = EndRefinement.begin() ; itr != EndRefinement.end() ; ++itr)
    {
        max = (int) itr->getNumIterations() > max ? itr->getNumIterations() : max ;
    }
    return max;
}


bool MeshRefinement::isRefinementIteration(const int &iteration){

    vector<AutoRef>::iterator itr;
    for (itr = this->AutoRefinement.begin() ; itr!= this->AutoRefinement.end() ; ++itr){
       // itr->printAutoref();
        if ( itr->isRefIter( iteration) ) return true;
    }
    return false;
}
