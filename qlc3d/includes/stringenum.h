#ifndef STRINGENUM_H
#define STRINGENUM_H
#include <globals.h>
#include <string>
#include <list>
#include <vector>
class StringEnum{
    /*!
      A class for simplifying the handling of string type settings.
      class is initialised using the constructor taking in key name and
      possible values.
      */
    std::string key;
    std::list< std::string> values;
    std::vector<int> int_values;
    bool parseValues(const std::string &values_in);
    StringEnum(){}
public:
    //StringEnum(const char* key, const char *val_in, const int* ints_in); //val_in MUST BE A COMMA DELIMITED LIST WITH NO SPACES
    StringEnum(const char *key_in, const char *val_in);
    bool containsValue(const std::string &val) const;
    int getValueIndex(const std::string &val) const;
    void printErrorMessage(const std::string &sin)const;
    void setIntValues(std::vector<int> & ivals);
};

#endif
