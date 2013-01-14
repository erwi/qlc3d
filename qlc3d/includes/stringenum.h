#ifndef STRINGENUM_H
#define STRINGENUM_H
#include <globals.h>
#include <string>
#include <list>
class StringEnum{
    /*!
      A class for simplifying the handling of string type settings.
      class is initialised using the constructor taking in key name and
      possible values.
      */
    std::string key;
    std::list< std::string> values;

    bool parseValues(const std::string &values_in);
    StringEnum(){}
public:
    StringEnum(const char* key, const char *val_in); //val_in MUST BE A COMMA DELIMITED LIST WITH NO SPACES
    bool containsValue(const std::string &val) const;
    int getValueIndex(const std::string &val) const;
    void printErrorMessage(const std::string &sin)const;
};

#endif
