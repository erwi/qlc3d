#include <box.h>
#include <stdio.h>
#include <stdlib.h>
//default constructor
Box::Box(){
	Type = "Normal";

	Params[0] = 1;
	Params[1] = 0;

	X[0] = 0; X[1] = 0;
	Y[0] = 0; Y[1] = 0;
	Z[0] = 0; Z[1] = 0;


	Tilt[0]     = 0; Tilt[1]    = 0;
	Twist[0]    = 0; Twist[1]   = 0;
}

void Box::printBox(){
	printf("\tType: ");
	printf("\t%s",Type.c_str());
	printf("\n\tparams : %f, %f\n",Params[0],Params[1]);
	printf("\tX\t= [%1.1f, %1.1f]\n\tY\t= [%1.1f, %1.1f]\n\tZ\t= [%1.1f, %1.1f]\n",X[0],X[1],Y[0],Y[1],Z[0],Z[1]);
	printf("\tTilt\t= [%1.1f, %1.1f]\n",Tilt[0],Tilt[1]);
	printf("\tTwist\t= [%1.1f, %1.1f]\n",Twist[0],Twist[1]);
}
void Box::setParams(std::vector<double> p){
if (p.size() == 2){
        Params[0] = p[0];
        Params[1] = p[1];
    }
    else{
    std::cout<< "error, Box::setParams, invalid Params length - bye!"<< std::endl;
    exit(1);
    }
}
void Box::setX(std:: vector<double> x ){
    if (x.size() == 2){
        X[0] = x[0];
        X[1] = x[1];
    }
    else{
    std::cout<< "error, Box::setX, invalid X length - bye!"<< std::endl;
    exit(1);
    }
}

void Box::setY(std:: vector<double> y ){
    if (y.size() == 2){
        Y[0] = y[0];
        Y[1] = y[1];
    }
    else{
    std::cout<< "error, Box::setY, invalid Y length - bye!"<< std::endl;
    exit(1);
    }
}

void Box::setZ(std:: vector<double> z ){
    if (z.size() == 2){
        Z[0] = z[0];
        Z[1] = z[1];
    }
    else{
    std::cout<< "error, Box::setZ, invalid Z length - bye!"<< std::endl;
    exit(1);
    }
}

void Box::setTilt(std:: vector<double> tlt ){
    if (tlt.size() == 2){

	Tilt[0] = tlt[0];
        Tilt[1] = tlt[1];
    }
    else{
    std::cout<< "error, Box::setTilt, invalid Tilt length - bye!"<< std::endl;
    exit(1);
    }
}
void Box::setTwist(std:: vector<double> twt ){
    if (twt.size() == 2){
        Twist[0] = twt[0];
        Twist[1] = twt[1];
    }
    else{
    std::cout<< "error, Box::setX, invalid Twist length - bye!"<< std::endl;
    exit(1);
    }
}
bool Box::isInBox(double *coords){

    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    //  test if coordinate is outside of box and return false
    if ( ( x < this->X[0] ) || ( x > this->X[1]) ) return false; // if smaller than minimum or larger than maximum...
    if ( ( y < this->Y[0] ) || ( y > this->Y[1]) ) return false;
    if ( ( z < this->Z[0] ) || ( z > this->Z[1]) ) return false;
    // otherwise return true
    return true;
}


Boxes::Boxes(){
		n_Boxes = 0;

}


Boxes::~Boxes(){

	std::vector<Box*>::iterator itr;

	for(itr = box.begin() ; itr!= box.end() ; itr++)
		delete (*itr);

}

void Boxes::addBox()
{
	box.push_back(new Box);
	n_Boxes++;

}

void Boxes::addBox(Box* b){
    box.push_back(b);
    n_Boxes ++;
}

void Boxes::printBoxes()
{
	std::vector<Box*>::iterator i;
	printf("%i Boxes:\n",n_Boxes);
	for (int i= 0 ; i < n_Boxes; i++)
	{
		printf("Box%i :\n",i+1);
		box[i]->printBox();

	}

}

void Boxes::WriteBoxes(FILE* fid)
{
	if (fid != NULL)
	{
		fprintf(fid, "#===============================\n");
		fprintf(fid, "#     INITIAL LC ORIENTATION\n");
		fprintf(fid, "#===============================\n\n");

		char str[50];
		for (int i = 0 ; i < n_Boxes ; i ++ )
		{
			sprintf(str,"BOX%i",i+1);


			fprintf(fid,"\t%s.Type = %s\n",str,   box[i]->Type.c_str() );
			fprintf(fid,"\t%s.Params = [%2.4f , %2.4f ]\n",str,    box[i]->Params[0] , box[i]->Params[1] );//Params
			fprintf(fid,"\t%s.X = [%2.4f, %2.4f]\n", str, box[i]->X[0] , box[i]->X[1]); //X
			fprintf(fid,"\t%s.Y = [%2.4f, %2.4f]\n", str, box[i]->Y[0] , box[i]->Y[1]); //Y
			fprintf(fid,"\t%s.Z = [%2.4f, %2.4f]\n", str, box[i]->Z[0] , box[i]->Z[1]); //Z
			fprintf(fid,"\t%s.Tilt  = [%2.4f, %2.4f]\n",str, box[i]->Tilt[0] , box[i]->Tilt[1]); //Tilt
			fprintf(fid,"\t%s.Twist = [%2.4f, %2.4f]\n\n", str, box[i]->Twist[0] , box[i]->Twist[1]);//Twist
		}



	}
}
