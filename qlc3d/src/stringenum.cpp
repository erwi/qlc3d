#include <stringenum.h>
#include <iostream>
#include <algorithm>

StringEnum::StringEnum(const char *key_in, const char* val_in)
{
    //val_in MUST BE A COMMA DELIMITED LIST WITH NO SPACES
    key = key_in;
    parseValues(val_in);
}

bool StringEnum::parseValues(const std::string &values_in)
{
    std::string line = values_in; // working copy

    while (line.size() )
    {
        for (idx i = 0 ; i < line.size() ; i++)
        {
            char ch = line.at(i);
            if (ch == ',')
            {
                std::string subs = line.substr(0,i);
                line.erase(0,i+1);
                values.push_back(subs);
                break;
            }

            if ((ch == '\0') || (i == line.size()-1))
            {
                std::string subs = line.substr(0,i+1);
                values.push_back(subs);

                line.erase(0,i+1); // REMOVE JUST ADDED SUBSTRING
                break;
            }


        }
    }
    return true;
}


int StringEnum::getValueIndex(const std::string &testval) const
{
    // LOOP OVER EACH VALUE AND COMPARE WITH VALID ONES
    // RETURNS INDEX TO FOUND VALUE, OR -1 IF NOT FOUND
    // ALL-LOWERCASE TEMPORARIES ARE USED, i.e. SEARCH IS
    // CASE INSENSITIVE


    std::string tval = testval; // working copy
    std::transform(tval.begin(), tval.end(), tval.begin(), ::tolower);

    std::list<std::string> ::const_iterator itr = values.begin();
    for (int retVal = 0; itr != values.end() ; itr++, retVal++)
    {
        std::string tval2 = *itr;
        std::transform(tval2.begin(), tval2.end(),tval2.begin(), ::tolower);

        if (tval2.compare(tval) == 0 )
        {
             return retVal;
        }
    }
    return -1; // Value not found, return negative index
}

bool StringEnum::containsValue(const std::string &testval) const
{
    std::string tval = testval; // working copy
    std::transform(tval.begin(), tval.end(), tval.begin(), ::tolower);

    std::list<std::string> ::const_iterator itr = values.begin();
    for (; itr != values.end() ; itr++)
    {
        std::string tval2 = *itr;
        std::transform(tval2.begin(), tval2.end(),tval2.begin(), ::tolower);

        if (tval2.compare(tval) == 0 ) // IF FOUND
        {
             return true;
        }
    }
    return false; // NOT FOUND
}


void StringEnum::printErrorMessage(const std::string &s_in)const
{
    std::cout << "Error specifying value for key : \"" << this->key <<"\""<< std::endl;
    std::cout << "Valid values are: " << std::endl;
    std::list<std::string> ::const_iterator itr = values.begin();
    for (; itr!=values.end() ; itr++)
    {
        std::cout <<"\t"<< *itr << std::endl;
    }

    std::cout << "Found \"" << s_in << "\" instead. Bye!" << std::endl;


}
