#include <refinfo.h>


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
    refIter_(0)

{
    setType(Type);
    //printRefInfo();
}

RefInfo::RefInfo(const RefInfo &other):
    type_(other.type_),
    iter_(other.iter_),
    time_(other.time_),
    refIter_(other.refIter_)
{

}



void RefInfo::setValues(std::vector<double> &values)
{
    values_.clear();
    values_.insert( values_.begin(), values.begin(), values.end() );
    this->setRefIter();
}

void RefInfo::printRefInfo(FILE *fid) const
{
    fprintf(fid, "Iteration = %li\n", this->iter_);
    fprintf(fid, "Time = %e\n", this->time_);
}
