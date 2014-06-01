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
    std::vector<std::string> values;
    std::vector<int> enumValues;   // enum values can be signed/unsigned

    StringEnum(){}
public:

    StringEnum(const std::string &key, const std::vector<std::string>& validStrings) {
        /*!
         * key is key name in settings file, used in preinting error message.
         * validStrings is array of valid string values
         * The values of the unum are assumed to be consecutive and starting at 0.
         */
        this->key = key;
        this->values = validStrings;
        // populate with consecutive values for enum
        for (int i = 0; i < (int) validStrings.size(); i++)
            enumValues.push_back(i);
    }

    StringEnum(const std::string &key,
               const std::vector<std::string>& validStrings,
               const std::vector<int>& enumValues) {
        /*!
         * key is key name in settings file, used in preinting error message.
         * validStrings is array of valid string values.
         * enum_values are enum values, useful in cases when not consecutive
         */
        // make sure enum names size matches number of custom enum numerical
        // values.
        if (validStrings.size() != enumValues.size()) {
            std::cerr << "error in " << __PRETTY_FUNCTION__ << std::endl;
        }

        this->key = key;
        this->values = validStrings;
        this->enumValues = enumValues;
    }


    EnumType getEnumValue(const std::string &sval) const {
   /*! Returns enum that corresponds to input string. If
    *  input string is not valid, an exception is thrown instead
    */
        //Make a working copy of input and make it all lower case
        std::string tval = sval; // working copy
        std::transform(tval.begin(), tval.end(), tval.begin(), ::tolower);
        // loop over all valid strings
        size_t enum_index = 0;
        for ( ; enum_index < values.size(); enum_index++){
            std::string valid_value = values[enum_index];
            std::transform(valid_value.begin(), valid_value.end(),valid_value.begin(), ::tolower);
            // If input value is valid, can cast index to enum and return early
            if (valid_value.compare(tval) == 0 ) {
                return  static_cast<EnumType> (enumValues[enum_index]);
            }
        }
        // input string did not belong to set of valid strings
        throw std::exception();
    }

    void printErrorMessage(const std::string &s_in) const {
        std::cerr << "Error specifying value for key : \"" << this->key <<"\""<< std::endl;
        std::cerr << "Valid options are: " << std::endl;
        for (auto &v : values)
            std::cerr << "\t" << v << std::endl;
        std::cerr << "Found \"" << s_in << "\" instead. Bye!" << std::endl;
    }
};

#endif
