#ifndef STRINGENUM_H
#define STRINGENUM_H
#include <globals.h>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdexcept>
template <typename EnumType>
class StringEnum{
    /*!
      A class for simplifying the handling of string type settings.
      class is initialised using the constructor taking in key name and
      possible values.
      */
    std::string key;
    std::list< std::string> values;
    std::vector<int> int_values;
    bool parseValues(const std::string &values_in){
        // BREAKS UP COMMA SEPARATED STRINGS LIST
        std::string line = values_in; // working copy
        while (line.size() ){
            for (idx i = 0 ; i < line.size() ; i++){
                char ch = line.at(i);
                if (ch == ','){
                    std::string subs = line.substr(0,i);
                    line.erase(0,i+1);
                    values.push_back(subs);
                    break;
                }

                if ((ch == '\0') || (i == line.size()-1)){
                    std::string subs = line.substr(0,i+1);
                    values.push_back(subs);
                    line.erase(0,i+1); // REMOVE JUST ADDED SUBSTRING
                    break;
                }
            }
        }
        return true;
    }
    StringEnum(){}
public:

    StringEnum(const char *key_in, const char *val_in, const int* iVals = NULL){
        // CONSTRUCTOR WITH DEFAULT INTEGER VALUE LIST
        // val_in MUST BE A COMMA DELIMITED LIST OF STRINGS WITH NO SPACES
        // iVals MUST BE EITHER NULL OR POINTER TO INTEGER ARRAY OF SAME LENGTH AS val_in
        key = key_in;
        parseValues(val_in);
        // IF DEFAULT, NO VARIABLE ARGUMENTS

        for (int i = 0 ; i < (int) this->values.size() ; i++){
            if (iVals==NULL){// IFF NULL, USE DEFAULTS
                int_values.push_back(i);
            }
            else {
                int_values.push_back(iVals[i]);
            }
        }
    }
    EnumType getEnumValue(const std::string &sval) const {
        // LOOP OVER EACH STRING VALUE AND COMPARE WITH VALID ONES
        // RETURNS NUMBER ENUMERATOR WHEN VALID STRING IS FOUND
        // OTHERWISE THROWS EXCEPTION
        std::string tval = sval; // working copy
        std::transform(tval.begin(), tval.end(), tval.begin(), ::tolower);

        std::list<std::string> ::const_iterator itr = values.begin();
        for (int retVal = 0; itr != values.end() ; itr++, retVal++){
            std::string tval2 = *itr;
            std::transform(tval2.begin(), tval2.end(),tval2.begin(), ::tolower);

            if (tval2.compare(tval) == 0 ){
                 return  ( static_cast<EnumType> (int_values[retVal]) );
            }
        }
        throw std::exception();
    }
    void printErrorMessage(const std::string &s_in)const{
        std::cout << "Error specifying value for key : \"" << this->key <<"\""<< std::endl;
        std::cout << "Valid options are: " << std::endl;
        std::list<std::string> ::const_iterator itr = values.begin();
        for (; itr!=values.end() ; itr++){
            std::cout <<"\t"<< *itr << std::endl;
        }
        std::cout << "Found \"" << s_in << "\" instead. Bye!" << std::endl;
    }

};

#endif
