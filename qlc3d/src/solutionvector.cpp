#include <solutionvector.h>
#include <material_numbers.h>

SolutionVector::~SolutionVector(){
    ClearFixed();

    if (Values){ free(Values); }
    Values = NULL;

    if (Elim){	free(Elim);}
    Elim = NULL;
    if (EquNodes!=NULL){ free(EquNodes); }
    EquNodes = NULL;
}
void SolutionVector::ClearFixed(){
    if ( FixedNodes && FixedValues){
        free(FixedNodes);
        free(FixedValues);
        FixedNodes = NULL;
        FixedValues = NULL;
        setnFixed(0);
    }
		
    if (IsFixed){
        free(IsFixed);
        IsFixed = NULL;
    }
}

SolutionVector& SolutionVector::operator=(const SolutionVector& r){
    if (this==&r)
	return *this;

    // CLEAR ALL DATA
    ClearAll();

// REALLOCATE AND COPY
    nDoF	= r.nDoF;
    nFixed	= r.nFixed;
    nDimensions = r.nDimensions;
    nFreeNodes  = r.nFreeNodes;
// VALUES
    if (r.Values){
        Values = (double*) malloc( nDoF * nDimensions * sizeof(double) );
        memcpy ( Values, r.Values , nDoF * nDimensions * sizeof(double) );
    }

// FIXED NODED
    if (r.FixedValues){
        FixedValues = (double*) malloc( nFixed * nDimensions * sizeof(double) );
        memcpy(FixedValues, r.FixedValues, nFixed*nDimensions * sizeof(double) );
    }
    if (r.FixedNodes){
        FixedNodes  = (int*)    malloc( nFixed * nDimensions * sizeof(int)    );
        memcpy(FixedNodes, r.FixedNodes, nFixed*nDimensions * sizeof(int) );
    }

// PERIODIC EQU NODES AND ELIM
    if (r.Elim){
        Elim 	= (int*)    malloc(nDoF * nDimensions *  sizeof(int));
        memcpy( Elim, r.Elim , nDoF*nDimensions * sizeof(int) );
    }
    if (r.EquNodes){
        EquNodes 	= (int*)    malloc(nDoF * nDimensions *  sizeof(int));
        memcpy( EquNodes, r.EquNodes, nDoF * nDimensions * sizeof (int) );
    }
    if (r.IsFixed){
        IsFixed	= (bool*)   malloc(nDoF * nDimensions *  sizeof(bool) );
        memcpy( IsFixed, r.IsFixed, nDoF * nDimensions * sizeof(bool) );
    }

    return *this;
}

SolutionVector::SolutionVector(){
    nDoF = 0;
    nFreeNodes = 0;
    nFixed = 0;
    IsVector = false;
    nDimensions = 0;
    FixedNodes 		= NULL;
    FixedValues 	= NULL;
    Values 			= NULL;
    IsFixed 		= NULL;
    Elim 			= NULL;
    EquNodes 		= NULL;
}
SolutionVector::SolutionVector(int np){

    IsVector 		= false;
    Elim 			= NULL;
    EquNodes 		= NULL;
    nDimensions 	= 1;
    nDoF 			= np;
    nFreeNodes 		= np;
    nFixed 			= 0;
    FixedNodes 		= NULL;
    FixedValues 	= NULL;
        if ((Values = (double*)malloc(np*sizeof(double)))== NULL){
            printf("error - SolutionVector::SolutionVector(int np) - could not allocate memory\n");
            exit(1);
        }
        else{
            setValuesTo((double) 0);
        }
        IsFixed = NULL;
	
}
SolutionVector::SolutionVector(int np, int dim)
{
	//SolutionVector();
        Elim                = NULL;
        EquNodes            = NULL;
        IsVector            = true;
        nDimensions         = dim;
        nDoF                = np;
        nFreeNodes          = np;
        nFixed              = 0;
        FixedNodes          = NULL;
        FixedValues         = NULL;
        IsFixed             = NULL;
        Values              = NULL;
	Allocate( (unsigned int) np, (unsigned int) dim );

}

void SolutionVector::Allocate(const unsigned int &np, const unsigned int &ndim){
        /*
    if(Values != NULL) free( Values );

	Values = NULL;

	Values = (double*) malloc ( ndim*np*sizeof(double) );
	setValuesTo( 0.0 );
	setnDoF( np );
	setnDimensions( ndim );
	nFreeNodes = np;
	if ( Values == NULL){
		printf("error - SolutionVector::Aloocate - could not do it - bye!\n");
		exit(1);
	}
        */
    Resize(np , ndim); // OOPS. TWO NAMES, ONE FUNCTION, SOUNDS FAMILIAR...
}

void SolutionVector::setnDoF(int n)		{	nDoF = n;}
void SolutionVector::setnFixed(int n)		{	nFixed = n;}
void SolutionVector::setnDimensions(int n)	{	nDimensions = n;}
void SolutionVector::setValuesTo(const double& value)
{// sets all Values to value
	#pragma omp parallel for
	for (int i = 0 ; i < nDoF ; i ++)
		Values[i] = value;
}// end void setValuesTo

void SolutionVector::setValuesTo(const double *values){
	unsigned int n = nDoF * nDimensions;
	#pragma omp parallel for
	for (unsigned int i = 0 ; i < n ; i++){
		this->Values[i] = values[i];

	}
}
void SolutionVector::setValuesTo(const SolutionVector &other){
    this->setValuesTo( other.Values );
}
void SolutionVector::Resize(const unsigned int &n, const unsigned int &dim){
    ClearAll();
    nDoF = n;
    nFreeNodes = n;     // all nodes are free until set fixed
    nDimensions = dim;
    Values = (double*) malloc(nDoF * nDimensions * sizeof(double) );
}

void SolutionVector::ClearAll(){
    // CLEAR ALL DATA
    nDoF = 0;
    nFixed = 0;
    nDimensions = 0;
    nFreeNodes = 0;

    if (IsFixed)    free(IsFixed);      IsFixed     = NULL;
    if (Elim)	    free(Elim);         Elim        = NULL;
    if (EquNodes)   free(EquNodes);     EquNodes    = NULL;
    if (FixedNodes) free(FixedNodes);   FixedNodes  = NULL;
    if (FixedValues)free(FixedValues);  FixedValues = NULL;
    if (Values)	    free(Values);       Values      = NULL;

}

void SolutionVector::setFixedNodesQ(Alignment* alignment, Mesh* e)
{	//! sets fixed nodes index list and value lists, assuming values are currently correct
    //! FixedNodes is index list to all nodes that are fixied
    //! FixedValues are the corresponding individual values that are fixed
 
// 1. First get index to all strong anchoring nodes 
	vector <int> ind_to_nodes;
	vector <int> temp_index;
	vector <int>::iterator itr;
	
	for (int i = 0 ; i < alignment->getnSurfaces() ; i ++){

        if (alignment->IsStrong(i)){

            e->FindIndexToMaterialNodes( (i+1) * MAT_FIXLC1 , &temp_index );
            for ( itr = temp_index.begin() ; itr != temp_index.end(); itr ++){
				ind_to_nodes.push_back(*itr);

                }
		}
	}// end for i

//2. allocate memory for index and fixed values	
	if (FixedNodes != NULL ) free(FixedNodes);
	if (FixedValues != NULL) free(FixedValues);
	
	FixedNodes = (int*) malloc(getnDimensions() * ind_to_nodes.size() * sizeof(int) );
		if (FixedNodes == NULL){
			printf("error - SolutionVector::setFixedNodesQ() - could not allocate memory for fixed nodes, bye!");
			exit(1);
		}
	FixedValues = (double*) malloc( getnDimensions() * ind_to_nodes.size() * sizeof(double) );
		if (FixedValues == NULL){
			printf("error - SolutionVector::setFixedNodesQ() - could not allocate memory for fixed values, bye!");
			exit(1);
		}
		
// 3. copy index and values to arrays
	int i = 0;
	nFixed = ind_to_nodes.size();

      //  nFreeNodes = nDoF - nFixed;
        for (itr = ind_to_nodes.begin() ; itr != ind_to_nodes.end() ; itr++ ){
		FixedNodes[i] = *itr;
		FixedNodes[i + 1*nFixed] = *itr + 1*getnDoF();
		FixedNodes[i + 2*nFixed] = *itr + 2*getnDoF();
		FixedNodes[i + 3*nFixed] = *itr + 3*getnDoF();
		FixedNodes[i + 4*nFixed] = *itr + 4*getnDoF();
		
		FixedValues[i]            = getValue(*itr, 0); //q1
		FixedValues[i + 1*nFixed] = getValue(*itr, 1); //q2
		FixedValues[i + 2*nFixed] = getValue(*itr, 2); //q3;
		FixedValues[i + 3*nFixed] = getValue(*itr, 3); //q4;
		FixedValues[i + 4*nFixed] = getValue(*itr, 4); //q5;
		i++;
	}
	

	setBooleanFixedNodeList();
}
void SolutionVector::setFixedNodes(vector<int> *Material, vector<double> *val ,int *Elem,int *Mat,int nElem,int nNodes)
{// sets fixed node lists wit node numbers and corresponding fixed values
		// *Material = fixed material numbers defined in settings file
		// *val = values corresponding to *Material
		// *Elem = array of element node numbers
		//*Mat  = array of material numbers corresponding to elem
		// nElem and nNodes = dimensions of element node arrays 
	
        if (Material->size()!=val->size()){
                printf("\nerror-solutionvector::SetFixedNodes\nnumber of fixed values must equal number of materials - bye!");
		exit(1);
	}
	
	nFixed = 0; // reset number of fixed nodes
	if (FixedNodes != NULL) free(FixedNodes);	// these may have been set earlier, and should be reset
	if (FixedValues!= NULL) free(FixedValues);
	
	vector<int> vec_FixedNodes;
	vector<double> vec_FixedValues;
	
	
	//int n_fix_mat = Material->size();
	//printf("%i fixed materials\n",n_fix_mat);
	
	vector<int>::iterator i;
	vector<int>::iterator j; 
	vector<double>::iterator Val; // fixed value
	
        for ( i = Material->begin(), Val =val->begin() ; i != Material->end() ; i++,Val++ ){

		// find index to all elements of material i
		vector<int> fix_elem; 
		vector<int> fix_nodes;
		for (int x = 0 ; x < nElem ; x ++)
			{
			//printf(" %i == %i\n",*i,Mat[x]);
				if ( (Mat[x] & 31*MAT_ELECTRODE1) == *i )  // bitwise testing of material number 
				{
				fix_elem.push_back(x);
				//printf("found\n");
				}
			}
		
		// ...then get all nodes from the fixed elements
		for (j = fix_elem.begin(); j != fix_elem.end() ; j ++)
			{
			for (int x = 0 ; x < nNodes ; x ++)
				{
				//printf("%i ",*j);
				fix_nodes.push_back(Elem[(*j)*nNodes+x]);
				}
			}
		// remove repeated nodes from list
		sort(fix_nodes.begin(), fix_nodes.end());
		vector<int>::iterator sorted;
		sorted = unique(fix_nodes.begin(), fix_nodes.end());
		
		//add nodes and values to final vector
		for (j = fix_nodes.begin() ; j != sorted ; j++)
			{
				//printf("%i ",*j);
				vec_FixedNodes.push_back(*j);
				vec_FixedValues.push_back(*Val);
			}
	}

		// Allocate memory for fixed node lists
		nFixed = vec_FixedNodes.size();
		FixedNodes = (int*)malloc(nFixed * sizeof(int));
		FixedValues= (double*)malloc(nFixed * sizeof(double));

               // nFreeNodes = nDoF - nFixed; // update free nodes count
		// Copy values form vectors to memory
		int c= 0;
		for (i = vec_FixedNodes.begin(), Val = vec_FixedValues.begin(), c = 0;
			i != vec_FixedNodes.end() ; i ++, Val++, c++)
		{
			FixedNodes[c] = *i;
			FixedValues[c] = *Val;
		}
	
	
	
}
void SolutionVector::setFixedNodesPot(Electrodes* electrodes, Mesh* surface_mesh)
{// sets fixed nodes for potentials
	
	if (nFixed>0)
	{

		if (FixedNodes!=NULL) {free(FixedNodes); FixedNodes = NULL;}
		if (FixedValues!=NULL) {free(FixedValues); FixedValues = NULL;}
		if (IsFixed!=NULL) {free(IsFixed); IsFixed = NULL;}
		nFixed = 0;
	}
	
	for (int i = 0 ; i < electrodes->nElectrodes ; i++){

		AddFixed((i+1)*MAT_ELECTRODE1, electrodes->E[i]->Potential[0], surface_mesh);
		
	}	
	
	//set boolean list of fixed nodes
	setBooleanFixedNodeList();
	
}

void SolutionVector::setFixedNodesPot(Electrodes* electrodes,
                                      Mesh* surface_mesh,
									  double CurrentTime){
// sets fixed nodes for potentials, taking into account voltage waveforms and current time
// this method is horrible and need to be rewritten
	if (nFixed>0){ //IF RE-SETTING, CLEAR OLD
		if (FixedNodes) {free(FixedNodes); FixedNodes = NULL;}
		if (FixedValues) {free(FixedValues); FixedValues = NULL;}
		if (IsFixed) {free(IsFixed); IsFixed = NULL;}
        nFixed = 0;
    }
    for (int i = 0 ; i < electrodes->nElectrodes ; i++){
        // find index to current time in voltage waveform
        int indx =0;
        for ( int j = 0 ; j < electrodes->E[i]->getnTimes() ; j ++){
			if ( fabs(electrodes->E[i]->Time[j] - CurrentTime) < 1e-15){
                indx = j;
            }
        }

        AddFixed( (i+1)*MAT_ELECTRODE1, electrodes->E[i]->Potential[indx], surface_mesh);
    }

    setBooleanFixedNodeList();
}
// end void setFixedNodesPot
void SolutionVector::setBooleanFixedNodeList(){

	if (IsFixed != NULL) { free(IsFixed);}
	
		IsFixed = (bool*) malloc(nDimensions * nDoF * sizeof(bool));
		memset(IsFixed , false , nDimensions * nDoF * sizeof(bool)); // set all to false
		
	for (int i = 0 ; i < getnFixed() *nDimensions ; i ++) // then set only fixed nodes to true
		IsFixed[FixedNodes[i]]=true;

}
	

void SolutionVector::PrintFixedNodes(){//debugging
    printf("number of fixed nodes is : %i\n",nFixed);
    int i;
    for (i=0; i<nFixed*getnDimensions(); i++)
        printf("fixed node %i, node number %i, value %1.3f\n",i,FixedNodes[i],FixedValues[i]);

}
void SolutionVector::PrintValues()
{
	for (int i = 0 ; i < nDoF *nDimensions; i++)
	{
		printf("%i = %f\n", i, Values[i]);
	}

}
void SolutionVector::PrintElim(){
    printf("\nPrinting for dimension 0 only!:\n");
    for (int i = 0 ; i < nDoF ; i++)
        printf("Elim[%i] = %i\n",i,Elim[i]);
}
void SolutionVector::PrintEquNodes(){
    printf("\nEquNodes:\n");
    for (int i = 0 ; i < nDoF ; i ++) {
        printf("EquNodes[%i] = %i \n",i, EquNodes[i]);
    }
}
void SolutionVector::setToFixedValues()
{
	if ( (FixedNodes == NULL) || (FixedValues == NULL) ){
		printf("error - SolutionVector::setToFixedValues, NULL pointer - bye!\n");
		exit(1);
	}
	
	int i;
	for (i = 0 ; i < nFixed * getnDimensions(); i ++){
		Values[FixedNodes[i]] = FixedValues[i ];
	}

}

void SolutionVector::setPeriodicEquNodes(Geometry* geom){
/*! Generates equivalent node indexes for periodic boundaries.
  */
	if (Elim != NULL) free(Elim);						// allocate memory for equivalent nodes
	if (EquNodes != NULL) free(EquNodes);
	Elim 		= (int*) malloc(nDoF*nDimensions*sizeof(int));	
	EquNodes 	= (int*) malloc(nDoF*nDimensions*sizeof(int));
	
        memset(Elim    , 0 , nDoF*nDimensions*sizeof(int) );
        memset(EquNodes, 0 , nDoF*nDimensions*sizeof(int) );

        for (int i = 0; i < nDoF*nDimensions ; i++ ){ // set to 1,2,3....
		Elim[i] = i;
               // EquNodes[i] = i;
        }
        if (!geom->getleft_right_is_periodic() && !geom->gettop_bottom_is_periodic() && !geom->getfront_back_is_periodic() ){
		return; // no periodic boundaries, can return
	}
// NEED TO CONSIDER 3 DIFFERENT CASES, DEPENDING ON TYPE OF PERIODICITY OF MODELLING WINDOW

	double eps = 1e-5; // accuracy for coordinate comparisons
	
	 
	

//CASE 1 FRONT BACK IS PERIODIC ONLY
    if (geom->getfront_back_is_periodic() && !geom->getleft_right_is_periodic() && !geom->gettop_bottom_is_periodic() ){
			
    //SEPARATE NODES INTO TWO LISTS FOR FRONT AND BACK SURFACES
        list <int> front;
        list <int> back;

        for (int i = 0 ; i < nDoF ; i ++ ){
            if ( geom->getAbsYDist( i , geom->getYmin() ) <= eps ) // check if node i is on front surface
                { front.push_back(i) ; }
                else
                if ( geom->getAbsYDist( i , geom->getYmax() ) <= eps ) // check if node i is on back surface
                    { back.push_back(i); }
        }//end for i

        /// MAKE SURE EQUAL NUMBER OF NODES HAVE BEEN FOUND ON BOTH SURFACES
        if (front.size() != back.size() ){
            printf("error - SolutionVector::setPeriosdicEquNodes\n");
            printf("front and back surfaces do not have same number of nodes\n");
            printf("front = %i , back = %i - bye!\n " , (int) front.size() , (int) back.size() );
            exit(1);
        }
	
	///SEARCH FOR NODE EQUIVALENCIES BY COMPARING X AND Y COORDINATES
	///BACK NODES MAP TO FRONT NODES		
        list <int>:: iterator B;
        list <int>:: iterator F;
        int fc, bc; // cdebug counters
        for (F = front.begin(),fc = 0; F!=front.end() ; F++, fc++){
            bool found = false;
            int ind_n  = 0;
            double dist = 1000000;
            for (B = back.begin() , bc = 0; B!= back.end() ; B++, bc ++){
			
                //compare x and z coordinates
                double xdist = geom->getAbsXDist(*F , geom->getpX(*B) );
                double zdist = geom->getAbsZDist(*F , geom->getpZ(*B) );
                double tdist = xdist*xdist + zdist*zdist;

                if (tdist < dist){ // debug info, keep track of nearest found node
                    dist = tdist; // nearest distance
                    ind_n = *B;  // index to nearest distance
                }
				
                if (tdist < eps*eps){ // compare squared distances
                    Elim[*B] = *F;
                    found = true;
                    break;
                }
            }// end for B
            if (!found)	{
                printf("error - matching front/back periodic nodes not found - bye!\n");
                printf("front = p[%i] = %f, %f, %f\n", *F, geom->getpX(*F) , geom->getpY(*F) , geom->getpZ(*F) );
                printf("nearest back = p[%i] = %f, %f, %f\n", ind_n , geom->getpX(ind_n) , geom->getpY(ind_n) , geom->getpZ(ind_n) );
                printf("distance = %f, fc = %i, bc = %i\n", dist, fc, bc);
                exit(1);
            }
        }//end for F

	}
	else
// CASE 2 FRONT-BACK AND LEFT-RIGHT ARE PERIODIC

        if (geom->getfront_back_is_periodic() && geom->getleft_right_is_periodic() && !geom->gettop_bottom_is_periodic() ){


		// separate nodes into 8 lists, 4 x corners left/right and front/back planes 
		list <int> corn0; //x = 0, y = 0
		list <int> corn1; //x = 0, y = max
		list <int> corn2; //x = max, y = max
		list <int> corn3; //x = max, y = 0
		list <int> front; //x = 0 
		list <int> back;  //x = max
		list <int> right; //y = max
		list <int> left;  //y = 0
		//printf("nDoF = %i\n",nDoF);
                for (int i = 0 ; i < nDoF ; i++){ // loop over all nodes and insert to correct list


			if ( ( geom->getAbsXDist( i , geom->getXmin() ) <= eps ) && ( geom->getAbsYDist( i , geom->getYmin() ) <= eps ) ) // corn0
				{corn0.push_back(i);}
			else
			if ( ( geom->getAbsXDist( i , geom->getXmin() ) <= eps ) && ( geom->getAbsYDist( i , geom->getYmax() ) <= eps ) ) // corn1
				{corn1.push_back(i);}
			else
			if ( ( geom->getAbsXDist( i , geom->getXmax() ) <= eps ) && ( geom->getAbsYDist( i , geom->getYmax() ) <= eps ) ) // corn2
				{corn2.push_back(i);}
			else
			if ( ( geom->getAbsXDist( i , geom->getXmax() ) <= eps ) && ( geom->getAbsYDist( i , geom->getYmin() ) <= eps ) ) // corn3
				{corn3.push_back(i);}
			else
			if ( geom->getAbsYDist( i , geom->getYmin() ) <= eps ) // front surface
				{front.push_back(i);}
			else
			if ( geom->getAbsYDist( i , geom->getYmax() ) <= eps ) // back surface
				{back.push_back(i);}
			else 
			if ( geom->getAbsXDist( i , geom->getXmin() ) <= eps ) // left surface
				{left.push_back(i);}
			else
			if ( geom->getAbsXDist( i , geom->getXmax() ) <= eps ) // right surface
				{right.push_back(i);}
		

		}

		if ( ( corn0.size() != corn1.size() ) || ( corn0.size() != corn2.size() ) || ( corn0.size() != corn3.size() )  || ( left.size() != right.size() ) || (front.size() != back.size() ) )
		{
					printf("error, different number of periodic boundary nodes\n");
					printf("corners 0,1,2,3 = %i,%i,%i,%i\n", (int) corn0.size() , (int) corn1.size() , (int) corn2.size() , (int) corn3.size() );
					printf("front/back , left/right = %i/%i , %i/%i \n", (int) front.size() , (int) back.size() , (int) left.size() , (int) right.size() );
					exit(1);
		}
		else // OK
		{
					// debug messages
					//printf(" %f, %f, %f, %f\n", geom->getXmin() , geom->getXmax() , geom->getYmin() , geom->getYmax() );
					//printf("OK np = %i, corn0.size() = %i\n", geom->getnp() , (int) corn0.size() );
					//printf(" left,rigth = %i , %i\n", (int) left.size() , (int) right.size() );
					//printf(" front,back = %i , %i\n", (int) front.size() , (int) back.size() );
					//printf("size of c0,c1,c2,c3 =  %i, %i, %i, %i,\n", (int) corn0.size(), (int) corn1.size() , (int) corn2.size() , (int) corn3.size() );
					//nFreeNodes = nDoF - 3*corn0.size() - left.size() - back.size() ;
					//printf("nFreeNodes = %i\n",nFreeNodes);
		}
		
	// START FILLING IN ELIM ARRAY WITH CORNER NODES
	// corn 0 nodes stay unchanged but corn1,2,3 map to them
		list<int> :: iterator c0;
		list<int> :: iterator c1;
		list<int> :: iterator c2;
		list<int> :: iterator c3;
		
                for (c0 = corn0.begin() ; c0 != corn0.end() ; c0++ ){ // outer corner loop - corn0
                        bool found = false;
			double dist = 0;
                        for (c1 = corn1.begin() ; c1 != corn1.end() ; c1++){ // inner corner loop - corn1
                                dist = geom->getAbsZDist( *c0 , geom->getpZ(*c1) );
				if ( dist <= eps ) // same z coordinate -> match!
				{
					Elim[*c1] = *c0;
					found = true;
					break;
				}
			}
			if (!found)
			{printf("error - corner node 1 not found - bye\n"); exit(1); }
			
			
			found = false;
                        for ( c2 = corn2.begin() ; c2 != corn2.end() ; c2++ ){ // inner corner loop - corn2
                                dist = geom->getAbsZDist( *c0 , geom->getpZ(*c2) );
				if ( dist <= eps ) // same z cordinate - > match!
				{
					Elim[*c2] = *c0;
					found = true;
					break;
				} 
			}
			if (!found)
			{printf("error - corner node 2 not found - bye\n"); exit(1); }
		
			found = false;
                        for ( c3 = corn3.begin() ; c3 != corn3.end() ; c3++ ){ // inner corner loop - corn2
                                dist = geom->getAbsZDist( *c0 , geom->getpZ(*c3) );
                                if ( dist <= eps ){ // same z cordinate - > match!
                                        Elim[*c3] = *c0;
					found = true;
					break;
				} 
			}
			if (!found)
			{printf("error - corner node 3 not found - bye\n"); exit(1); }
				
		}

	// THEN FILL IN LEFT/RIGHT NODES
	// LEFT IS UNCHANGED AND RIGHT SIDE MAPS TO LEFT
		list<int> :: iterator L;
		list<int> :: iterator R;
		int lc, rc; // left and right list counters - mainly for debugging
                for ( L = left.begin(), lc = 0 ; L!= left.end() ; L++ , lc++){ // outer left/right side loop

			bool found = false;
			int ind_n  = 0;
			double dist = 1000000;
                        for (R = right.begin() , rc = 0; R!= right.end() ; R++ , rc++){
                            // compare y and z coordinates
                            double ydist = geom->getAbsYDist( *L , geom->getpY( *R ) );
                            double zdist = geom->getAbsZDist( *L , geom->getpZ( *R ) );
                            double tdist = ydist*ydist + zdist*zdist;
				
                            if (tdist < dist){// debud info, keep track of nearest found node
                                dist = tdist; // nearest distance
                                ind_n = *R;  // index to nearest distance
                            }

                            if ( tdist <= eps*eps ){ // compares distances squared
                                Elim[*R] = *L;
                                found = true;
                                break;
                            }
			}
                        if(!found){
				printf("error - matching left/right periodic nodes not found - bye\n"); 
				printf("left  = p[%i] = %f, %f, %f\n", *L, geom->getpX(*L) , geom->getpY(*L) , geom->getpZ(*L) );
				printf("nearest right = p[%i] = %f, %f, %f\n", ind_n, geom->getpX(ind_n) , geom->getpY(ind_n) , geom->getpZ(ind_n) );
				printf("distance = %f, lc = %i , rc = %i\n", dist, lc, rc);
				exit(1);
                        }
		
		}// end for L
		//printf("all left/right nodes found\n");
		
	// FINALLY FILL IN FRONT/BACK NODES
	// FRONT IS UNCHANGED, BACK MAPS TO FRONT
		list <int>:: iterator B;
		list <int>:: iterator F;
		int fc, bc; // cdebug counters 
		for (F = front.begin(),fc = 0; F!=front.end() ; F++, fc++)
		{
			bool found = false;
			int ind_n  = 0;
			double dist = 1000000;
			for (B = back.begin() , bc = 0; B!= back.end() ; B++, bc ++)
			{
			
				//compare x and z coordinates
				double xdist = geom->getAbsXDist(*F , geom->getpX(*B) );
				double zdist = geom->getAbsZDist(*F , geom->getpZ(*B) );
				double tdist = xdist*xdist + zdist*zdist;
				if (tdist < dist) // debud info, keep track of nearest found node
					{
						dist = tdist; // nearest distance
						ind_n = *B;  // index to nearest distance
					}
				
				
				if (tdist < eps*eps) // compare squared distances
				{
					Elim[*B] = *F;
					found = true;
					break;
				}
			}// end for B
			if (!found)
				{
					printf("error - matching front/back periodic nodes not found - bye!\n");
					printf("front = p[%i] = %f, %f, %f\n", *F, geom->getpX(*F) , geom->getpY(*F) , geom->getpZ(*F) );
					printf("nearest back = p[%i] = %f, %f, %f\n", ind_n , geom->getpX(ind_n) , geom->getpY(ind_n) , geom->getpZ(ind_n) );
					printf("distance = %f, fc = %i, bc = %i\n", dist, fc, bc);
					exit(1);
				}
		
		}//end for F
		
			
		//exit(0);
	}
// CASE 3 FRONT-BACK, LEFT-RIGHT AND TOP-BOTTOM ARE PERIODIC
	// not implemented yet
	else
        if (geom->getfront_back_is_periodic() && geom->getleft_right_is_periodic() && geom->gettop_bottom_is_periodic() ){
		printf("front/back, left/right and top/bottom are periodic, but this has not been implemented yet - bye!");
		exit(1);
	}
	
	
	
	
	
        for (int ii = 0; ii < nDoF ; ii ++){
                        EquNodes[ii] = Elim[ii];
				//printf("EquNodes[%i] = %i \n",ii+1, EquNodes[ii]+1);
        }

        nFreeNodes = 0;
		
	// adjust indices
	// find "wrong nodes" (nodes that are equivalent only)
	vector <int> wrong;
	vector <int> wrongval;
	for (int i = 0 ; i < nDoF ; i++)
	{	
		if ( Elim[i]!= i )	
			{
			wrong.push_back(i);
			wrongval.push_back(Elim[i]);
			}
	}	
	
	// reduce all indices after an "equivalent" node  	
	for (unsigned int i = 0 ; i < wrong.size() ; i++)
		{
		for (int j = wrong[i]; j <nDoF; j++ )
			Elim[j] --;
		}
	
	// insert back all "equivalent" indices  (although they are wrong at this point)
	for (unsigned int i = 0 ; i < wrong.size() ; i ++)
		Elim[wrong[i]] = wrongval[i];


 
		
	// fix all "equivalent" indices to point to correct nodes		
        for (unsigned int i = 0 ; i < wrong.size() ; i ++){
		Elim[wrong[i]] = Elim[Elim[wrong[i]]];
	}
	
	nFreeNodes = 0;
	for (int i = 0 ; i < nDoF ; i ++)
		if (nFreeNodes<Elim[i]) nFreeNodes = Elim[i];
	
	nFreeNodes ++;
	
	if (nDimensions > 1) // if more thand one dimension, copy index to other dimensions
	{
		for (int j = 1; j < nDimensions ; j ++)
		{
			for (int i = 0 ; i < nDoF ; i ++ )
				{
				Elim[j*nDoF + i] = Elim[i] + j*nFreeNodes;
				
				}
			
		}
	}


}
void SolutionVector::setValue(const unsigned int& n, const unsigned int& dim, const double& val)
{
	#ifdef DEBUG
	if (( (int) n >= getnDoF() ) || ( (int) dim >= getnDimensions() ) )
	{
		printf("error - SolutionVector::setValue(n, dim, val) - n = %i and dim = %i, bye!\n",n,dim);
		printf(", when\n");
		printf(" this->nDimensions = %i\n", this->nDimensions );
		printf(" this->nDoF = %i\n", this->nDoF );
		exit(1);
	}
	#endif
	Values[n + dim*getnDoF()] = val;
	
}

void SolutionVector::AddFixed(int mat, double val, Mesh *mesh){
    vector <int> ind_p; // index to all nodes of material mat
	

    mesh->FindIndexToMaterialNodes(mat,&ind_p);

	
// Allocate memory for old + new fixed nodes and values	
	int NewFixedSize = ind_p.size() + nFixed; // number of old + new fixed nodes

	int* NewFixedNodes = (int*) malloc(NewFixedSize*sizeof(int)); // allocate enough memory for old + new
	if (NewFixedNodes == NULL){
		printf("error - SolutionVector::AddFixed(int,double,Mesh*) - could not allocate memory for fixed nodes, bye!");
		exit(1);
	}
	
	double* NewFixedValues = (double*) malloc(NewFixedSize*sizeof(double));
	if (NewFixedValues == NULL)	{
		printf("error - SolutionVector::AddFixed(int,double,Mesh*) - could not allocate memry for fixed values, bye!");
		exit(1);
	}

// Copy old fixed nodes and values and add new ones 	
	vector <int>::iterator itr;
	itr = ind_p.begin();
	for (int i = 0 ; i < NewFixedSize ; i++){
		if (i<nFixed){ // copy old ones
			NewFixedNodes[i] = FixedNodes[i];
			NewFixedValues[i]= FixedValues[i];
		}
		else{ // add new ones
			NewFixedNodes[i] = *itr; //  = iterator to index to material nodes
			itr++;
		
			NewFixedValues[i] = val;
		}
	}// end for i
	
// free old arrays and redirect pointer to New arrays
	if (FixedValues!=NULL) free(FixedValues);
	if (FixedNodes !=NULL) free(FixedNodes);
	
	FixedValues = NewFixedValues;
	FixedNodes  = NewFixedNodes;
	
	NewFixedValues = NULL;
	NewFixedNodes  = NULL; // necessary to keep only single pointer to memory location?
	
	setnFixed(NewFixedSize); // update counter

	
}
void SolutionVector::EnforceEquNodes()
{
	if (EquNodes==NULL)
	{
		printf("error - SolutionVector::EnforceEquNodes(), Equnodes = NULL - bye\n");
		exit(1);
	}
	for (int i =0 ; i < nDoF ; i ++)
		for (int d = 0 ; d < 5 ; d++)
		Values[EquNodes[i]+d*nDoF] = Values[i+d*nDoF];


}









