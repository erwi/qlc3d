#include <refinfo.h>
#include <assert.h>

void RefInfo::setType(std::string s)
{
// TESTS WHETHER "REFINEMENT.Type" STRING IN SETTNGS FILE IS VALID.
// IF YES, SETS type_ ENUMERATOR FOR THIS OBJECT
// IF NOT, TERMINATES PROGRAM

    // MAKE LOWERCASE
    std::transform(s.begin(), s.end(), s.begin(), std::ptr_fun<int, int>(std::tolower));

    if ( ! s.compare("change") )
    {
        type_ = Change;
        return;
    }
    else
        if ( ! s.compare("sphere") )
        {
            type_ = Sphere;
            return;
        }
    else
    {
        std::cout<< "error - bad refinement type: "<< s << " bye!" << std::endl;
        exit(1);
    }
}
void RefInfo::setRefIter()
{
// DETERMINES NUMBER OF REFINEMENT ITERATIONS THIS OBJECT DESCRIBES
// AND SETS THE refIter_ VARIABLE

// DEPENDING ON THE type_, NUMBER OF REFINEMENT ITERATIONS IS DETERMINED BY
// DIFFERENT CONDITIONS
    switch(type_)
    {
    case(Change):
    {
        refIter_= (int) values_.size();
        break;
    }
    case (Sphere):
    {
        refIter_ = (int) values_.size();
        // MAKE SURE NOT EMPTY AND EQUAL NUMBER OF COORDS HAVE BEEN SPECIFIED
        break;
    }
    default:
        printf("error in %s, unhandled refinement type - bye\n", __func__);
        exit(1);
    }
}


RefInfo::RefInfo(const std::string& Type):
    type_(None),
    iter_(0),
    time_(0),
    refIter_(0),
    X_(1,0),
    Y_(1,0),
    Z_(1,0)

{
    setType(Type);

}

RefInfo::RefInfo(const RefInfo &other):
    type_(other.type_),
    iter_(other.iter_),
    time_(other.time_),
    refIter_(other.refIter_),
    values_(other.values_),
    X_(other.X_),
    Y_(other.Y_),
    Z_(other.Z_)
{

}



void RefInfo::setValues(std::vector<double> &values)
{
    values_.clear();
    values_.insert( values_.begin(), values.begin(), values.end() );
    this->setRefIter();
}


void RefInfo::setCoords(const std::vector<double> &x,
                        const std::vector<double> &y,
                        const std::vector<double> &z)
{
    X_.clear();
    Y_.clear();
    Z_.clear();

    X_.insert(X_.begin(), x.begin(), x.end() );
    Y_.insert(Y_.begin(), y.begin(), y.end() );
    Z_.insert(Z_.begin(), z.begin(), z.end() );

    this->setRefIter();

}

void RefInfo::getCoords(std::vector<double> &x,
                        std::vector<double> &y,
                        std::vector<double> &z) const
{
    x = X_;
    y = Y_;
    z = Z_;
}

void RefInfo::getCoord(double &x, double &y, double &z) const
{
    assert( X_.size() );
    assert( Y_.size() );
    assert( Z_.size() );

    x = X_[0];
    y = Y_[0];
    z = Z_[0];
}

double RefInfo::getValue(const size_t i) const
{
#ifdef DEBUG
    assert ( i <= values_.size() );
#endif
    return values_[i];
}

void RefInfo::printRefInfo(FILE *fid) const
{
    fprintf(fid, "Iteration = %li\n", this->iter_);
    fprintf(fid, "Time = %e\n", this->time_);
    fprintf(fid, "Type = %i\n", this->type_);

    for (size_t i = 0 ; i < values_.size() ; i++)
        fprintf(fid,"values_[%u] = %e\n", i , values_[i]);

}

void RefInfo::validate(const RefInfo &refinfo)
{
#define ERRORMSG printf("error, bad or missing REFINEMENT data - bye!\n")
    switch ( refinfo.getType() )
    {
    case( RefInfo::Change ):
        // MAKE SURE VALUES EXIST
        if (!refinfo.values_.size() )
        {
            ERRORMSG;
            exit(1);
        }
        break;
    case( RefInfo::Sphere ):
        // MAKE SURE VALUES AND COORDINATES EXIST
        if ( ( !refinfo.values_.size() ) ||
             ( !refinfo.X_.size() ) ||
             ( !refinfo.Y_.size() ) ||
             ( !refinfo.Z_.size() ) )
        {
             ERRORMSG;
             exit(1);
        }

        break;
    default:
        printf("error in %s, unknonwn refinement type - bye!\n");
        exit(1);
    }


}
