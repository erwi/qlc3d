#include <solutionvector.h>
#include <material_numbers.h>

const double SolutionVector::BIGNUM = 1e99;

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
    if ( FixedNodes && FixedValues)
    {
        free(FixedNodes);
        free(FixedValues);
        FixedNodes = NULL;
        FixedValues = NULL;
        setnFixed(0);
    }
		
    if (IsFixed)
    {
        free(IsFixed);
        IsFixed = NULL;
    }
    if (FixedNodeMaterial)
    {
        free(FixedNodeMaterial);
        FixedNodeMaterial = NULL;
    }
}

SolutionVector& SolutionVector::operator=(const SolutionVector& r)
{
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
    //ClearFixed();
    if (r.FixedValues)
    {
        FixedValues = (double*) malloc( nFixed * nDimensions * sizeof(double) );
        memcpy(FixedValues, r.FixedValues, nFixed*nDimensions * sizeof(double) );
    }

    if (r.FixedNodes)
    {
        FixedNodes  = (int*)    malloc( nFixed * nDimensions * sizeof(int)    );
        memcpy(FixedNodes, r.FixedNodes, nFixed*nDimensions * sizeof(int) );
    }

    if (r.FixedNodeMaterial)
    {
        FixedNodeMaterial = (int*) malloc( nFixed*sizeof(int) );
        memcpy( FixedNodeMaterial, r.FixedNodeMaterial, nFixed*sizeof(int) );
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

SolutionVector::SolutionVector():
    nDoF(0),
    nFixed(0),
    nDimensions(0),
    FixedNodeMaterial(NULL),
    IsFixed(NULL),
    Elim(NULL),
    EquNodes(NULL),
    nFreeNodes(0),
    IsVector(false),  // <-- is this used anywhere ??
    FixedNodes(NULL),
    FixedValues(NULL),
    Values(NULL)
{

}
SolutionVector::SolutionVector(int np):
    nDoF(np),
    nFixed(0),
    nDimensions(1),
    FixedNodeMaterial(NULL),
    IsFixed(NULL),
    Elim(NULL),
    EquNodes(NULL),
    nFreeNodes(np),
    IsVector(false),  // <-- is this used anywhere ??
    FixedNodes(NULL),
    FixedValues(NULL),
    Values(NULL)

{
    Allocate( (unsigned int) np , (unsigned int) nDimensions );
}
SolutionVector::SolutionVector(int np, int dim):
    nDoF(np),
    nFixed(0),
    nDimensions(dim),
    FixedNodeMaterial(NULL),
    IsFixed(NULL),
    Elim(NULL),
    EquNodes(NULL),
    nFreeNodes(np),
    IsVector(true),  // <-- is this used anywhere ??
    FixedNodes(NULL),
    FixedValues(NULL),
    Values(NULL)
{
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
{
// ALL VALUES ARE SET TO THAT OF THE INPUT VALUE
    #pragma omp parallel for
    for (size_t i = 0 ; i < (size_t) nDoF ; i ++)
        Values[i] = value;
}// end void setValuesTo

void SolutionVector::setValuesTo(const double *values)
{
// ALL VALUES ARE SET TO THOSE OF THE INPUT VECTOR
// VECTOR LENGHTS MUCH MATCH, NO CHECKING IS PERFORMED HERE!!
    size_t n = nDoF * nDimensions;
    #pragma omp parallel for
    for (size_t i = 0 ; i < n ; i++)
        this->Values[i] = values[i];
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

void SolutionVector::ClearAll()
{
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
    if (FixedNodeMaterial) free(FixedNodeMaterial); FixedNodeMaterial = NULL;
}

void SolutionVector::setFixedNodesQ(Alignment* alignment, Mesh* e)
{	/*! sets fixed nodes index list and value lists, assuming values are currently correct
		FixedNodes is index list to all nodes that are fixied
		FixedValues are the corresponding individual values that are fixed
	*/
 
// 1. First get index to all strong anchoring nodes 



    vector <int> ind_to_nodes;
    for (int i = 0 ; i < alignment->getnSurfaces() ; i ++)
    {
        if (alignment->IsStrong(i))
        {
            printf("FIXLC%i is strong\n", i+1 );
            vector <int> temp_index;
            e->FindIndexToMaterialNodes( (i+1) * MAT_FIXLC1 , &temp_index );
            //for ( itr = temp_index.begin() ; itr != temp_index.end(); itr ++)
            //{
            //    ind_to_nodes.push_back(*itr);
            //}

            ind_to_nodes.insert(ind_to_nodes.end(), temp_index.begin(), temp_index.end() );
        }
        else
        {
            printf("FIXLC%i is not strong\n", i+1);
            printf("%s\n", alignment->surface[i]->getAnchoringType().c_str() );
        }
    }// end for i


    // INDEX TO FIXED NODES MAY CONTAIN DUPLICATED ENTRIES, IF TWO FIXED
    // LC SURFACES ARE NEX TO EACHOTHER. REMOVE THESE

    sort(ind_to_nodes.begin() , ind_to_nodes.end() );
    vector <int>::iterator itr;
    itr = unique(ind_to_nodes.begin() , ind_to_nodes.end() );
    ind_to_nodes.erase(itr, ind_to_nodes.end() );

//2. allocate memory for index and fixed values	
	if (FixedNodes != NULL ) free(FixedNodes);
	if (FixedValues != NULL) free(FixedValues);
        nFixed = ind_to_nodes.size();

        // IF NO FIXED NODES EXIST, LEAVE
        if (!nFixed )
        {
            setBooleanFixedNodeList();  // SET ALL TO NON-FIXED
            return;
        }


        FixedNodes = (int*) malloc(getnDimensions() * nFixed * sizeof(int) );
        if (FixedNodes == NULL){
			printf("error - SolutionVector::setFixedNodesQ() - could not allocate memory for fixed nodes, bye!");
			exit(1);
		}
        FixedValues = (double*) malloc( getnDimensions() * nFixed * sizeof(double) );
		if (FixedValues == NULL){
			printf("error - SolutionVector::setFixedNodesQ() - could not allocate memory for fixed values, bye!");
			exit(1);
		}
		
// 3. copy index and values to arrays. This assumes that current Q-tensor values are correct
		// and will be fixed at these values (frozen to current)
	int i = 0;
	nFixed = ind_to_nodes.size();

      //  nFreeNodes = nDoF - nFixed;
        for (itr = ind_to_nodes.begin() ; itr != ind_to_nodes.end() ; itr++ )
        {
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

void SolutionVector::allocateFixedNodesArrays(Geometry &geom)
{
// ALLOCATES ARRAYS FOR MANAGING FIXED NODES.
//
// CALCULATES NUMBER OF FIXED NODES nFixed
// SETS FIXED NODES INDEX ARRAY FixedNodes VALUES
// ALLOCATES MEMORY FOR FIXED VALUES, BUT DOES NOT SET VALUES
// SETS FIXED NODES BOOLEAN FLAG ARRAY (CALLS setBoleanFixedNodesList)




    // THIS IS AN INITIALIZATION ONLY FUNCTION.
    // RESET ALL FIRST, IF RESIZING SOLUTION VECTOR
    // (E.G. AFTER MESH REFIENEMENT)
    if (nFixed>0)
    {
        printf("error in %s, fixed arrays are already initialised - bye!\n",__func__);
        exit(1);
    }

    // SEPARATE INDEX TO ALL FIXED NODES
    // (SOME WILL BE REPEATED)
    Mesh &e = *geom.e;   // PTR TO SURFACE MESH
    std::vector<SolutionVectorNameSpace::node> fixed_nodes;
    for (int i = 0 ; i < e.getnElements() ; i++)
    {
        int mat = e.getMaterialNumber(i);
        size_t indE = MATNUM_TO_ELECTRODE_NUMBER((size_t) mat);

        if ( indE ) // IF ELECTRODE ELEMENT
        {
            for (int j = 0 ; j < e.getnNodes() ; j++)
                fixed_nodes.push_back(
                            SolutionVectorNameSpace::
                            node( e.getNode(i,j), mat ) );
        }
    }
    // REMOVE REPEATED NODE INDEXES
    sort(fixed_nodes.begin(), fixed_nodes.end() );
    std::vector<SolutionVectorNameSpace::node> :: iterator itr;
    itr = unique(fixed_nodes.begin(), fixed_nodes.end() );
    fixed_nodes.erase(itr, fixed_nodes.end() );

    nFixed = fixed_nodes.size();

    // CREATE FIXED NODES INDEXES AND VALUES ARRAYS
    FixedNodes = (int*) malloc( nFixed*sizeof(int) * nDimensions );
    FixedNodeMaterial = (int*) malloc( nFixed*sizeof(int) );
    FixedValues= (double*) malloc( nFixed*nDimensions*sizeof(double) );

    if ( (!FixedNodes) || (!FixedValues) || (!FixedNodeMaterial) )
    {
        printf("error in %s, malloc returned NULL - bye!\n", __func__);
        exit(1);
    }

    // SET FIXED NODES INDEX ARRAY VALUES
    for (int i = 0 ; i < nFixed ; i++)
    {
        for (int j = 0 ; j <nDimensions ; j++)
        {
            FixedNodes[i+j*nFixed] = fixed_nodes[i].nodenum;
        }
        FixedNodeMaterial[i] = fixed_nodes[i].mat;
    }


    // ALL FIXED VALUES TO 0
    memset(FixedValues, 0 , nFixed*nDimensions*sizeof(double) );

    setBooleanFixedNodeList();  // SET BOOLEAN FLAGS
}

void SolutionVector::setFixedNodesPot(Electrodes* electrodes)
{
    // SETS VALUES IN FixedValues

    for (int i = 0 ; i < nFixed ; i++)
    {
        int mat = FixedNodeMaterial[i];
        size_t indE = MATNUM_TO_ELECTRODE_NUMBER((size_t) mat );

        if (indE)
        {
            double pot = electrodes->getCurrentElectrodePotential( indE - 1 );
            FixedValues[i] = pot;
        }
    }

	
}

void SolutionVector::setFixedNodesPot(  Electrodes& electrodes,
                                        Mesh* surface_mesh,
                                        double CurrentTime)
{
// sets fixed nodes for potentials, taking into account voltage waveforms and current time
// this method is horrible and need to be rewritten
    if (nFixed>0)
    { //IF RE-SETTING, CLEAR OLD
        if (FixedNodes) {free(FixedNodes); FixedNodes = NULL;}
        if (FixedValues) {free(FixedValues); FixedValues = NULL;}
        if (IsFixed) {free(IsFixed); IsFixed = NULL;}
        nFixed = 0;
    }
    
	
    for (size_t i = 0 ; i < electrodes.getnElectrodes() ; i++) // for each electrode
    {
    // find index to current time in voltage waveform
        //int indx =-1; // start search with invalid index
        
        //for ( int j = 0 ; j < electrodes->E[i]->getnTimes() ; j ++) // for each switching time
        //{
        //    if ( fabs(electrodes->E[i]->Time[j] - CurrentTime) < 1e-15) // if switchin occuring now
        //    {
        //        indx = j;
        //        break;
        //    }
        //}

        // If a switching event was found, apply new boundary conditions and remove event
        //if (indx >= 0 )
        //{
        //    double pot = electrodes->E[i]->Potential[indx]; // value of new potential
        //    electrodes->E[i]->setCurrentPotential( pot );   //
        //}
        double pot = electrodes.getCurrentElectrodePotential( i );

        // set all fixed nodes for electrode i to its current potential
        AddFixed( (i+1)*MAT_ELECTRODE1,
                  pot ,
                  surface_mesh);

    }// end for each electrode

    setBooleanFixedNodeList();
}
// end void setFixedNodesPot


void SolutionVector::setBooleanFixedNodeList(){
// A BOOLEAN FLAG FOR EACH NODE, WHETHER IT IS FIXED OR NOT.
// MAYBE A PARAMETERS VARIABLE WIT BIT MASKS WOULD BE MORE EFFICIENT
// IF MANY "FLAG" ARRAYS ARE NEEDED
    if (IsFixed != NULL) { free(IsFixed); }
	
// ALLOCATE MEMORY FOR ARRAY
    size_t size = nDimensions * nDoF * sizeof( bool );
    IsFixed = (bool*) malloc( size );
    memset(IsFixed , false , size ); // set all to false

// SET VALUE TO TRUE/FALSE FOR EACH NODE
    for (int i = 0 ; i < getnFixed() *nDimensions ; i ++) // then set only fixed nodes to true
        IsFixed[FixedNodes[i]]=true;

}
	

void SolutionVector::PrintFixedNodes()
{
    //debugging
    printf("number of fixed nodes is : %i\n",nFixed);

    for ( int i = 0 ; i < nFixed ; i++ )
    {
        for (int j = 0 ; j < nDimensions ; j++)
        {
            printf("fixed node %i, node number %i, value %1.3f, ",i,FixedNodes[i],FixedValues[i]);
        }
        if(FixedNodeMaterial)
            printf("material %i", FixedNodeMaterial[i]);
        printf("\n");
    }
}
void SolutionVector::PrintValues()
{
	for (int i = 0 ; i < nDoF *nDimensions; i++)
	{
                printf("Value[%i] = %e\n", i, Values[i]);
	}

}
void SolutionVector::PrintElim(){

    if (!Elim)
    {
        printf("Elim array is NULL\n");
        return;
    }

    printf("\nPrinting for dimension 0 only!:\n");
    for (int i = 0 ; i < nDoF ; i++)
        printf("Elim[%i] = %i\n",i,Elim[i]);
}
void SolutionVector::PrintEquNodes(){
    if (!EquNodes)
    {
        printf("EquNodes array is NULL\n");
        return;
    }

    printf("\nEquNodes:\n");
    for (int i = 0 ; i < nDoF ; i ++) {
        printf("EquNodes[%i] = %i \n",i, EquNodes[i]);
    }
}
void SolutionVector::PrintIsFixed()
{
    printf("%i fixed, %i free nodes:\n", this->getnFixed(), this->getnFreeNodes() );
    for (int i = 0 ; i < this->getnDoF()*this->getnDimensions() ; i++)
    {
        printf("node[%i] = ", i);
        if (getIsFixed(i))
            printf("TRUE\n");
        else
            printf("FALSE\n");
    }

}

void SolutionVector::setToFixedValues()
{
// SETS ALL VALUES TO CORRECT FIXED VALUES

    // MAKE SURE ARRAYS HAVE BEEN INITIALISED
    if ( ( (FixedNodes == NULL) || (FixedValues == NULL)) &&
                ( nFixed > 0 )    )
    {
        printf("error - SolutionVector::setToFixedValues, NULL pointer - bye!\n");
        exit(1);
    }

    for (int i = 0 ; i < nFixed * getnDimensions(); i ++)
    {
        int ind = FixedNodes[i];
        double val = FixedValues[i];
        Values[ ind ] = val;
    }
}




void SolutionVector::setPeriodicEquNodes(Geometry* geom){
/*!
    SET VALUES IN THE ELIM ARRAY. THE ELIM ARRAY CONTAINS
    MAPPINGS FORM A NODE NUMBER TO ITS ACTUAL
    DEGREE OF FREEDOM (ITS POSITION IN THE GLOBAL MATRIX)
*/
    // IF NO PERIODIC NODES PRESENT, DON'T GENERATE EQUIVALENT NODES INDEXES
    if (!geom->getleft_right_is_periodic() &&
        !geom->gettop_bottom_is_periodic() &&
        !geom->getfront_back_is_periodic() &&
           ( this->nFixed == 0 )
            )
    {
        return; // no periodic boundaries, can return
    }


    if (Elim != NULL) free(Elim);						// allocate memory for equivalent nodes
    if (EquNodes != NULL) free(EquNodes);
    Elim        = (int*) malloc(nDoF*nDimensions*sizeof(int));

    for (int i = 0; i < nDoF*nDimensions ; i++ )
        Elim[i] = i; // set to 0,1,2,3....

	
// CREATE LIST OF PERIODIC NODE INDEXES
    vector <unsigned int> periNodes;
    geom->e->listNodesOfMaterial( periNodes, MAT_PERIODIC );


// NEED TO CONSIDER 3 DIFFERENT CASES, DEPENDING ON TYPE OF PERIODICITY OF MODELLING WINDOW
    double eps = 1e-5; // accuracy for coordinate comparisons
    double xmin = geom->getXmin(); // convenience shortcuts to geometry min and max dimensions
    double xmax = geom->getXmax();
    double ymin = geom->getYmin();
    double ymax = geom->getYmax();
    double zmin = geom->getZmin();
    double zmax = geom->getZmax();
    /// PROBABLY EVIL, BUT SO CONVENIENT...
    #define LEFT    ( geom->getAbsXDist(n, xmin) <= eps )
    #define RIGHT   ( geom->getAbsXDist(n, xmax) <= eps )
    #define FRONT   ( geom->getAbsYDist(n, ymin) <= eps )
    #define BACK    ( geom->getAbsYDist(n, ymax) <= eps )
    #define BOTTOM  ( geom->getAbsZDist(n, zmin) <= eps )
    #define TOP     ( geom->getAbsZDist(n, zmax) <= eps )


//CASE 1 FRONT BACK IS PERIODIC ONLY
    if (    geom->getfront_back_is_periodic() &&
            !geom->getleft_right_is_periodic() &&
            !geom->gettop_bottom_is_periodic() )
    {
			
    //SEPARATE NODES INTO TWO LISTS FOR FRONT AND BACK SURFACES
        list <int> front;
        list <int> back;

        for (size_t i = 0 ; i < periNodes.size() ; i ++ )
        {
            unsigned int n = periNodes[i];

            if ( geom->getAbsYDist( n , geom->getYmin() ) <= eps ) // check if node i is on front surface
                { front.push_back(n) ; }
            else
            if ( geom->getAbsYDist( n , geom->getYmax() ) <= eps ) // check if node i is on back surface
                { back.push_back(n); }

            else // ERROR
            {
                printf("error in %s, CASE 1 - bye \n", __func__ );
                exit(1);
            }
        }//end for i


        /// MAKE SURE EQUAL NUMBER OF NODES HAVE BEEN FOUND ON BOTH SURFACES
        if (front.size() != back.size() ){
            printf("error - SolutionVector::setPeriosdicEquNodes\n");
            printf("front and back surfaces do not have same number of nodes\n");
            printf("front = %i , back = %i - bye!\n " , (int) front.size() , (int) back.size() );
            exit(1);
        }
	
        // SEARCH FOR NODE EQUIVALENCIES BY COMPARING X AND Y COORDINATES
        // BACK NODES MAP TO FRONT NODES
        setFaceElim( front, back, Elim, 1, geom->getPtrTop() );
    }
	else
// CASE 2 FRONT-BACK AND LEFT-RIGHT ARE PERIODIC
        if (    geom->getfront_back_is_periodic() &&
                geom->getleft_right_is_periodic() &&
                !geom->gettop_bottom_is_periodic() )
        {

            // separate nodes into 8 lists, 4 x corners left/right and front/back planes
            list <int> corn0; //x = 0, y = 0
            list <int> corn1; //x = 0, y = max
            list <int> corn2; //x = max, y = max
            list <int> corn3; //x = max, y = 0
            list <int> front; //x = 0
            list <int> back;  //x = max
            list <int> right; //y = max
            list <int> left;  //y = 0

        for (size_t i = 0 ; i < periNodes.size() ; i++) // loop over all nodes and insert to correct list
        {
            int n = periNodes[i];
            if ( LEFT && FRONT ) // corn0
                {corn0.push_back(n);}
            else
            if ( LEFT && BACK ) // corn1
                {corn1.push_back(n);}
            else
            if ( RIGHT  &&  BACK ) // corn2
                {corn2.push_back(n);}
            else
            if ( RIGHT && FRONT ) // corn3
                {corn3.push_back(n);}
            else
            if ( FRONT ) // front surface
                {front.push_back(n);}
            else
            if ( BACK ) // back surface
                {back.push_back(n);}
            else
            if ( LEFT ) // left surface
                {left.push_back(n);}
            else
            if (RIGHT ) // right surface
                {right.push_back(n);}
        }

        if ( ( corn0.size() != corn1.size() ) ||
            ( corn0.size() != corn2.size() )  ||
            ( corn0.size() != corn3.size() )  ||
            ( left.size() != right.size() ) ||
            (front.size() != back.size() ) )
        {
            printf("error, different number of periodic boundary nodes\n");
            printf("corners 0,1,2,3 = %i,%i,%i,%i\n", (int) corn0.size() , (int) corn1.size() , (int) corn2.size() , (int) corn3.size() );
            printf("front/back , left/right = %i/%i , %i/%i \n", (int) front.size() , (int) back.size() , (int) left.size() , (int) right.size() );
            exit(1);
        }
        setCornerElim( corn0, corn1, corn2, corn3, Elim, 2, geom->getPtrTop() ); // vertical corners
        setFaceElim( left, right, Elim, 0, geom->getPtrTop() ); // left/right faces
        setFaceElim( front, back, Elim, 1, geom->getPtrTop() ); // front/back faces

    }// END CASE 2

// CASE 3 FRONT-BACK, LEFT-RIGHT AND TOP-BOTTOM ARE PERIODIC

	else
    if (geom->getfront_back_is_periodic() && 
        geom->getleft_right_is_periodic() &&
        geom->gettop_bottom_is_periodic() )
	{
            // separate nodes into lists, 12 x edges left/right, front/back and top/bottom planes
            // Vertical corners along Z

            // Additionally, 7 corner nodes must point to bottom left (origin xmin,ymin,zmin) corner

            list <int> corn0; //x = 0, y = 0
            list <int> corn1; //x = 0, y = max
            list <int> corn2; //x = max, y = max
            list <int> corn3; //x = max, y = 0

            // Horizontal corners along X
            list <int> corna; // y = 0, z = 0
            list <int> cornb; // y = max, z = 0
            list <int> cornc; // y = max, z = max
            list <int> cornd; // y = 0, z = max

            // Horizontal corners along Y
            list <int> cornA; // x = 0, z = 0
            list <int> cornB; // x = max, z = 0
            list <int> cornC; // x = max, z = max
            list <int> cornD; // x = 0, z = max

            list <int> front; //x = 0
            list <int> back;  //x = max
            list <int> right; //y = max
            list <int> left;  //y = 0
            list <int> top;   //z = max
            list <int> bottom;//z = 0;



            // LOOP OVER ALL NODES AND INSERT TO CORRECT LIST


            int corner_nodes[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
            for (size_t i = 0 ; i < periNodes.size() ; i++)
            {
                int n = periNodes[i];

                // CORNER NODES TAKE PRECEDENCE OVER OTHER NODES
                // FRONT LEFT BOTTOM
                if ( FRONT && LEFT && BOTTOM )
                    corner_nodes[0] = n;
                else
                // FRONT RIGHT BOTTOM
                if ( FRONT && RIGHT && BOTTOM )
                    corner_nodes[1] = n;
                else
                // FRONT LEFT TOP
                if ( FRONT && LEFT && TOP )
                    corner_nodes[2] = n;
                else
                // FRONT RIGHT TOP
                if (FRONT && RIGHT && TOP)
                    corner_nodes[3] = n;
                else
                // BACK LEFT BOTTOM
                if (BACK && LEFT && BOTTOM)
                    corner_nodes[4] = n;
                else
                // BACK RIGHT BOTTOM
                if (BACK && RIGHT && BOTTOM)
                    corner_nodes[5] = n;
                else
                // BACK LEFT TOP
                if (BACK && LEFT && TOP)
                    corner_nodes[6] = n;
                else
                // BACK RIGHT TOP
                if (BACK && RIGHT && TOP)
                    corner_nodes[7] = n;



                else
                // EDGE NODES
                // 4 x Vertical Corners
                if ( LEFT && FRONT ) // corn0
                {corn0.push_back(n);}
                else
                if ( LEFT  && BACK ) // corn1
                {corn1.push_back(n);}
                else
                if (RIGHT  && BACK ) // corn2
                {corn2.push_back(n);}
                else
                if ( RIGHT && FRONT ) // corn3
                {corn3.push_back(n);}
                else

                // 4 x Horizontal along X
                if ( FRONT && BOTTOM)	// ymin and zmin
                    {corna.push_back(n);}
                else
                if (BACK && BOTTOM) // ymax and zmin
                    {cornb.push_back(n);}
                else
                if (BACK && TOP) // ymax and zmax
                    {cornc.push_back(n);}
                else
                if (FRONT && TOP)  // ymin and zmax
                    {cornd.push_back(n);}

                // 4 x Horizontal along Y

                else
                if( LEFT && BOTTOM)
                    {cornA.push_back(n);} // xmin and zmin
                else
                if( RIGHT && BOTTOM)
                    {cornB.push_back(n);} // xmax and zmin
                else
                if( RIGHT && TOP )
                    {cornC.push_back(n);} // xmax and zmax
                else
                if( LEFT && TOP)
                    {cornD.push_back(n);} // xmin and zmax

                else
                // FRONT/BACK, LEFT/RIGHT, TOP/BOTTOM FACES
                if ( FRONT ) // front surface
                    {front.push_back(n);}
                else
                if ( BACK ) // back surface
                    {back.push_back(n);}
                else
                if ( LEFT ) // left surface
                    {left.push_back(n);}
                else
                if ( RIGHT ) // right surface
                    {right.push_back(n);}
                else
                if ( BOTTOM ) // bottom surface
                    {bottom.push_back(n); }
                else
                if ( TOP ) // top surface
                    {top.push_back(n); }
		
            }// end for i, loop over all nodes

            // CHECK THAT OPPOSITE FACES HAVE EQUAL NUMBER OF NODES
            {// start dummy scope

                int mincorner = *min_element( corner_nodes, corner_nodes+8);
                if (mincorner < 0)
                {
                    printf(" error - corner nodes not found - bye!\n");
                    printf("indexes are = [%i,%i,%i,%i,%i,%i,%i,%i]\n", corner_nodes[0],corner_nodes[1],corner_nodes[2],corner_nodes[3],
                           corner_nodes[4],corner_nodes[5],corner_nodes[6],corner_nodes[7]);
                    exit(1);

                }


                if (top.size() != bottom.size() )
		{
                    printf("error - top and bottom surfaces do not match - bye!\n");
                    printf("sizes are top,bottom = %i,%i\n", top.size(), bottom.size() );
                    exit(1);
		}
		if (left.size() != right.size() )
		{
                    printf("error - left and right surfaces do not match - bye!\n");
                    exit(1);
		}
		if (front.size() != back.size() )
		{
                    printf("error - front and back surfaces do not match - bye!\n");
                    exit(1);
		}
		// CHECK ALL CORNERS HAVE CORRECT NUMBER OF NODES
		size_t s0, s1, s2, s3;
		s0 = corn0.size(); s1 = corn1.size(); s2 = corn2.size() ; s3 = corn3.size();
		if ( (s1!=s0) || (s2 != s0) || (s3!=s0) )
		{
                    printf("error - vertical corner node counts do not match\n");
                    exit(1);
		}
		s0 = corna.size(); s1 = cornb.size(); s2 = cornc.size(); s3 = cornd.size();
		if ( (s1!=s0) || (s2 != s0) || (s3!=s0) )
		{
                    printf("error - horizontal corner (along x) node counts do not match\n");
                    exit(1);
		}
		s0 = cornA.size(); s1 = cornB.size(); s2 = cornC.size(); s3 = cornD.size();
		if ( (s1!=s0) || (s2 != s0) || (s3!=s0) )
		{
                    printf("error - horizontal corner (along y) node counts do not match\n");
                    exit(1);
		}		
            }// end dummy scope


            // set corner nodes
            Elim[corner_nodes[1] ] = corner_nodes[0];
            Elim[corner_nodes[2] ] = corner_nodes[0];
            Elim[corner_nodes[3] ] = corner_nodes[0];
            Elim[corner_nodes[4] ] = corner_nodes[0];
            Elim[corner_nodes[5] ] = corner_nodes[0];
            Elim[corner_nodes[6] ] = corner_nodes[0];
            Elim[corner_nodes[7] ] = corner_nodes[0];



            // match edge nodes
            setCornerElim( corn0, corn1, corn2, corn3, Elim, 2, geom->getPtrTop() ); // vertical corners
            setCornerElim( corna, cornb, cornc, cornd, Elim, 0, geom->getPtrTop() ); // horiz. along X
            setCornerElim( cornA, cornB, cornC, cornD, Elim, 1, geom->getPtrTop() ); // horiz. along Y
            // faces nodes
            setFaceElim( left, right, Elim, 0, geom->getPtrTop() ); // left/right faces
            setFaceElim( front, back, Elim, 1, geom->getPtrTop() ); // front/back faces
            setFaceElim( bottom,top , Elim, 2, geom->getPtrTop() ); // top/bottom faces

            //printf("front/back, left/right and top/bottom are periodic, but this has not been implemented yet - bye!");
            //exit(1);
	}// end if 3 different periodicity cases


    // NODAL EQUIVALENCIES HAVE BEEN SET.
    // REPLACE DEPENDENT NODES WITH THEIR
    // INDEPENDENT EQUIVALENT NODES


    // MARK FIXED NODES. THSE WILL BE REMOVED FROM
    // FREE DEGREES OF FREEDOM
    for (int i = 0 ; i < nDoF ; i++ )
    {
        if ( this->getIsFixed(i) )
            Elim[i] = FIXED_NODE;
    }



    nFreeNodes = 0;

    std::vector <int> elim(nDoF, 0 );   // convenience copy of Elim
    elim.insert(elim.begin(), Elim, Elim+nDoF);
    elim.resize( nDoF );

    std::vector <int> elima(nDoF,0);    // Elim altered
    for (int i = 0 ; i < nDoF ; i++)    // SET TO 1,2,3...
        elima[i] = i;
    elima.resize(nDoF);

    // LOOP OVER EACH NODE. DECREASE INDEX TO ALL INDEPENDENT DOFs
    // THAT COME AFTER A DEPENDENT NODE (EQUIVALENT TO SHIFTING LEFT
    // ROWS/COLUMNS OF A MATRIX AFTER A COLUMN IS REMOVED)
    for (int i = 0 ; i < nDoF ; i++)
    {
        if (elim[i] != i)   // IF i'th NODE IS DEPENDENT
        {
            for (int j = i ; j < nDoF ; j++) // SHIFT DOWN ALL DOF INDEXES AFTER IT
            {
                elima[j] --;
            }
        }
    }

    // SET DEPENDENT VARAIBLE INDEXES TO POINT TO CORRECT
    // INDEPENDENT DOF
    for (int i = 0 ; i < nDoF ; i++) // SET CORRECT VALUES
    {
        if ( (elim[i]!=i) && (elim[i]!= FIXED_NODE) ) // IF i'th NODE IS DEPENDENT ( AND NOT FIXED)
        {
            elima[i] = elima[ elim[i] ]; // IT WILL DEPEND ON THE CORRECTED DOF INDEX
        }
        else
        if ( elim[i] == FIXED_NODE )    // KEEP FIXED NODE FLAGS
        {
            elima[i] = FIXED_NODE;
        }

    }

    // TOTAL NUMBER OF FREE DOFs THAT NEED TO BE SOLVED (PER DIMENSION)
    nFreeNodes = *max_element(elima.begin(), elima.end() ) + 1;

    // COPY BACK VALUES TO Elim ARRAY
    for (int i = 0 ; i < nDoF ; i++)
        Elim[i] = elima[i];

    // EXPAND Elim IF MORE THAN ONE DIMENSIONS PER NODE
    if (nDimensions > 1)
    {
        for (int j = 1; j < nDimensions ; j ++)
        {
            for (int i = 0 ; i < nDoF ; i ++ )
            {
                if (Elim[i] == FIXED_NODE)
                    Elim[j*nDoF+i] = FIXED_NODE;
                else
                    Elim[j*nDoF + i] = Elim[i] + j*nFreeNodes;
            }
        }
    }

}// end setPeriondicEquNodes

void SolutionVector::setCornerElim(
                        list <int>& corn0, // sets periodic equivalent nodes for 4 corners
                        list <int>& corn1, // corn1[i] = corn0[i]
                        list <int>& corn2, // corn2[i] = corn0[i]
                        list <int>& corn3, // corn3[i] = corn0[i]
                        int* Elim,
                        const int& dim, // direction of corner 0,1,2 -> x,y,z
                        double* p // coordinates
                        )
					 
{
/*! This sets Elim vector for periodic corner nodes. 
 * Nodes are arranged so that corn1[i] = corn2[i] = corn2[i] = corn0[i]
 * AbsDist is a function pointer to a function that calculates absolute distance along 
 * X, Y or Z axes (use geom->getAbsNDist, where N is X, Y, or Z).
 * getCoord is a function pointer to getting X,Y or Z coordinate (replaces geom->getpN)
 * */	
	
    list<int> :: iterator c0;
    list<int> :: iterator c1;
    list<int> :: iterator c2;
    list<int> :: iterator c3;
    double eps = 1e-5;

    for (c0 = corn0.begin() ; c0 != corn0.end() ; c0++ )
    { // outer corner loop - corn0
        bool found = false;
        double dist = 0;
            
        double p0 = p[(*c0) * 3 + dim ]; // dim coordinate of node c0
        int iNearest = -1;
        double minDist = BIGNUM;

        // COMPARISON WITH C1
        for (c1 = corn1.begin() ; c1 != corn1.end() ; c1++)
        { // inner corner loop - corn1
            double p1 = p[(*c1) * 3 + dim];
            dist = fabs( p0 - p1 );

            if ( dist <= minDist )
            {
                minDist = dist;
                iNearest = *c1;
            }


            if ( dist <= eps ) // if same coordinate -> match!
            {
                Elim[*c1] = *c0;
                found = true;
                break;
            }
        }
        if (!found)
        {

            printf("error - corner node 1 not found - bye\n");
            printf("corner node 0 = %i,\tat [%e,%e,%e]\n", *c0, p[*c0*3 + 0], p[*c0*3 + 1], p[*c0*3 + 2] );
            printf("nearest node  = %i,\tat [%e,%e,%e]\n", iNearest, p[iNearest*3 + 0], p[iNearest*3 + 1], p[iNearest*3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
			
	// COMPARISON WITH C2		
        found = false;
        iNearest = -1;
        minDist = BIGNUM;

        for ( c2 = corn2.begin() ; c2 != corn2.end() ; c2++ )
        { // inner corner loop - corn2
            double p2 = p[(*c2)*3 + dim ];
            dist = fabs( p0 - p2 );

            if ( dist <= minDist )
            {
                minDist = dist;
                iNearest = *c2;
            }

            if ( dist <= eps ) // if same cordinate -> match!
            {
                Elim[*c2] = *c0;
                found = true;
                break;
            }
        }

        if (!found)
        {
            printf("error - corner node 2 not found - bye\n");
            printf("corner node 0 = %i,\tat [%e,%e,%e]\n", *c0, p[*c0*3 + 0], p[*c0*3 + 1], p[*c0*3 + 2] );
            printf("nearest node  = %i,\tat [%e,%e,%e]\n", iNearest, p[iNearest*3 + 0], p[iNearest*3 + 1], p[iNearest*3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
	// COMPARISON WITH C3
        found = false;
        iNearest = -1;
        minDist = BIGNUM;

        for ( c3 = corn3.begin() ; c3 != corn3.end() ; c3++ )
        { // inner corner loop - corn2
            double p3 = p[(*c3)*3 + dim ];
            dist = fabs(p0 - p3 );

            if ( dist <= minDist )
            {
                minDist = dist;
                iNearest = *c1;
            }
            if ( dist <= eps )  // if same cordinate -> match!
            {
                Elim[*c3] = *c0;
                found = true;
                break;
            }
        }
        if (!found)
        {
            printf("error - corner node 3 not found - bye\n");
            printf("corner node 0 = %i,\tat [%e,%e,%e]\n", *c0, p[*c0*3 + 0], p[*c0*3 + 1], p[*c0*3 + 2] );
            printf("nearest node  = %i,\tat [%e,%e,%e]\n", iNearest, p[iNearest*3 + 0], p[iNearest*3 + 1], p[iNearest*3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
    }// end for loop over corn0 nodes c0
}// end void setCornerElim
void SolutionVector::setFaceElim( list <int>& face0, // face1[i] = face0[i]
				list <int>& face1,
				int* Elim,
				const int& norm, // face normal, 0,1,2 -> x,y,z
				double* p) // pointer to node coordinates
{
    // norm is face normal, need coordinates to vectors parallel to it
    int ind1 = ( norm + 1 ) % 3; // ind is a pre-calculated offeset to coordinate comparison in p.
    int ind2 = ( norm + 2 ) % 3; // e.g. norm = 0 -> ind1 = 1, ind2 = 2
    double eps = 1e-5; // accuracy of coordinate comparison
	
	    
    // SEARCH FOR NODE EQUIVALENCIES BY COMPARING COORDINATES
    // THAT ARE PERPENDICULAR TO FACE NORMAL
    list <int>:: iterator F0;    // FACE 0 NODES
    list <int>:: iterator F1;    // FACE 1 NODES
    int fc, bc; // debug counters
    for (F0 = face0.begin(), fc = 0; F0!=face0.end() ; F0++, fc++ ) // LOOP OVER FACE 0
    {
        bool found = false;
        int ind_n  = 0;     // index to neares (debug)
        double dist = 1000000;
        double f1 = p[3*(*F0) + ind1 ]; // coordinates of node F2 in face0
        double f2 = p[3*(*F0) + ind2 ];
			
        for (F1 = face1.begin() , bc = 0; F1!= face1.end() ; F1++, bc ++)
        {
            double fa = p[3*(*F1) + ind1]; // coordinates of node F1 in face 1
            double fb = p[3*(*F1) + ind2];
				
            //compare coordinates
            double dist1 = fabs(f1 - fa); // distances in plane
            double dist2 = fabs(f2 - fb);
            double tdist = dist1*dist1 + dist2*dist2;

            if (tdist < dist) // debug info only, keep track of nearest found node
            {
                dist = tdist; // nearest distance
                ind_n = *F1;  // index to nearest distance
            }
				
            if (tdist < eps*eps) // compare squared distances
            {
                Elim[*F1] = *F0;
                found = true;
                break;
            }
        }// end for B
        if (!found)
        {
            printf("error - matching periodic faces with normal %i - bye!\n", norm);
                //printf("front = p[%i] = %f, %f, %f\n", *F, geom->getpX(*F) , geom->getpY(*F) , geom->getpZ(*F) );
                //printf("nearest back = p[%i] = %f, %f, %f\n", ind_n , geom->getpX(ind_n) , geom->getpY(ind_n) , geom->getpZ(ind_n) );
                //printf("distance = %f, fc = %i, bc = %i\n", dist, fc, bc);
                exit(1);
        }
    }//end for F
//*/
}				

void SolutionVector::setValue(const unsigned int& n, const unsigned int& dim, const double& val)
{
	//#define DEBUG
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

void SolutionVector::AddFixed(int mat, double val, Mesh *mesh)
{
    vector <int> ind_p; // index to all nodes of material mat

    mesh->FindIndexToMaterialNodes(mat,&ind_p);


// Allocate memory for old + new fixed nodes and values	
    int NewFixedSize = ind_p.size() + nFixed; // number of old + new fixed nodes

    int* NewFixedNodes = (int*) malloc(NewFixedSize*sizeof(int)); // allocate enough memory for old + new
    if (NewFixedNodes == NULL)
    {
        printf("error - SolutionVector::AddFixed(int,double,Mesh*) - could not allocate memory for fixed nodes, bye!");
        exit(1);
    }
	
    double* NewFixedValues = (double*) malloc(NewFixedSize*sizeof(double));
    if (NewFixedValues == NULL)
    {
        printf("error - SolutionVector::AddFixed(int,double,Mesh*) - could not allocate memry for fixed values, bye!");
        exit(1);
    }

// Copy old fixed nodes and values and add new ones 	
    vector <int>::iterator itr;
    itr = ind_p.begin();
    for (int i = 0 ; i < NewFixedSize ; i++)
    {
        if (i<nFixed) // copy old ones
        {
            NewFixedNodes[i] = FixedNodes[i];
            NewFixedValues[i]= FixedValues[i];
        }
        else// add new ones
        {
                NewFixedNodes[i] = *itr; //  = iterator to index to material nodes
                itr++;
		
                NewFixedValues[i] = val;
        }
    }// end for i
	
// free old arrays and set pointer to New arrays
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
// MAKES SURE THAT VALUES ON PERIODIC SURFACES ARE OK.
// THIS MAY BE NEEDED e.g. AT THE START OF A SIMULATION
// OR TO AVOID ACCUMULATION OF NUMERICAL NOISE(?)

    // IF NO REORDERING OF DEGREES OF FREEDOM, CAN LEAVE
    if (nFreeNodes == nDoF)
    {
            return;
    }

    if (Elim == NULL)
    {
        printf("error - SolutionVector::EnforceEquNodes(), Elim = NULL - bye\n");
        exit(1);
    }

    // CALCULATES PERIODIC EQUIVALEN NODE FROM ELIM IN CASES WHERE
    // WHERE MORE THAN 1 DEGREE OF FREEDOM EXISTS (i.e. Q-TENSOR)
    for (int i =0 ; i < nDimensions ; i ++)
    {
        for (int j = 0 ; j < nDoF; j++)
        {
            int equDof = getEquNode(j);
            if ( equDof != FIXED_NODE )
            {
                int dep = i*nDoF + j;          // DEPENDENT NODE
                int indep = i*nDoF + equDof;   // EQUIVALENT INDEPENDENT NODE

                Values[ dep ] = Values[ indep ];
            }// END IF NOT FIXED NODE
        }
     }
}









