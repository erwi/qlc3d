#ifndef READER_H
#define READER_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <locale>
#define READER_SUCCESS      1
#define READER_NOT_FOUND    -2   // variable not found
#define READER_BAD_VALUE    -3   //
#define READER_BAD_FORMAT   -4   // missing "=" etc.
#define READER_BAD_FILE	    -5   // file not open etc.
using std::cout;
using std::endl;
class Reader
{
    public:
	std::fstream file;

	bool isIgnoreCase;
	bool allowSpaces;

        Reader();
        virtual ~Reader();
        bool openFile(const std::string& filename);
        void closeFile();
        bool gotoLine(const unsigned int& l); // uggly seek line number in txt files

        int readString(std::string var, std::string& val);

        int readDlmLine(std::vector <double>& val, std::string dlm  );

        int findStringLine(std::string str,		// finds line consisting of string str,
			   unsigned int offset = 0);	// returns line number with respect to offset

        int readNextAssignment(std::string& var,
                               std::string& val);// reads next variable = value assignment;

        int findVariableLine(const std::string& var,   // returns line next number where variable var is assigned
                             unsigned int offset = 0); // if offset ==0 , line number is w.r.t. beginning of file

        std::string getErrorString(int e);

        //
        //	FILE WRITING FUNCTIONS
        //
        int writeLine(const std::string& line, const int &ln = -1 ); // writes a line, if ln >= 0 line is written on that line


    template <class T>
    int readNumberArray(std::string var, std::vector<T>& val)
    {
        val.clear();
        file.seekg(0);
        if ( ! file.good() )
            return READER_BAD_FILE;

        int ln = this->findVariableLine(var); // FIND VARIABLE WE WANT TO READ
        if (ln>=0)
        {// IF FOUND
            this->gotoLine( (unsigned int) ln);
            std::string line;
            getline(file , line );
            line = this->extractRHS(line);

            std::vector <double> vec;
            int ret = this->stringToDoubleVec(line , vec ); // CONVERT STRING TO DOUBLE VECTOR

            // CONVERT DOUBLE VECTOR TO DESIRED FORMAT
            if (ret == READER_SUCCESS)
            {
                std::vector<double>::iterator itr;
                for (itr = vec.begin() ; itr != vec.end() ; itr++ )
                {
                    val.push_back((T) *itr );
                }
            }
            file.seekg(0);
            file.clear();
            return ret;
        } // END IF FOUND
        return READER_NOT_FOUND;
    }

        template <class T>
        int readNumber(std::string var, T& val){
	    this->file.seekg(0);

			if ( !file.good() )
				return READER_BAD_FILE;

			int ln = this->findVariableLine(var);

			if (ln >= 0){ // if variable definition found
				std::string line;
				this->gotoLine(ln);
				getline(file, line);
				line = this->extractRHS(line);

				double temp_dbl;
				if ( std::stringstream( line ) >> temp_dbl ){
					val = (T) temp_dbl;
					return READER_SUCCESS;
				}
				else
					return READER_BAD_VALUE; // variable found, but not numeric (?)
			}
			return READER_NOT_FOUND;
        }
        // end int readValue


   //template <class T>
   int readStringArray(std::string var, std::vector<std::string>& val);






    private:
        bool removeBlanks(std::string& str);// removes spaces
        bool removeBlanksAtBeginning(std::string& str);// removes spaces at beginning of line only
        bool removeBlanksAtEnd(std::string& str);// removes empty spaces at end of string
        bool removeComments(std::string& str);// removes anythign after a "#" character
        void stringToLowercase(std::string& str); // converts string to all lowercase
	std::string extractRHS(std::string line); // extracts Right Hands Side from a variable definition
	std::string extractLHS(std::string line);
        int stringToDoubleVec(std::string& str , std::vector <double>& vec);
        // SPLITS COMMA DELIMITED STRING TO VECTOR OF STRINGS
        int stringToStringVec(std::string& str , std::vector <std::string>& ret_str);
};




#endif // READER_H
