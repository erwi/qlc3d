#include <filesysfun.h>

bool setCurrentDirectory(const std::string& destdir)
{
// Changes current working directory to destdir.
// returns true/false if/not succesful.
#ifdef Linux
    if ( chdir( destdir.c_str() ) == 0 )
        return true;
    else
        return false;
#endif
}

std::string getCurrentDirectory()
{
#ifdef Linux
    char* buff = get_current_dir_name(); // allocates space to buffer
    std::string cwd = buff;
    free( buff );   // must free buffer

    return cwd;
#endif
}

bool dirExists(const std::string& dir)
{
    // checks whether directory dir exists, returning true/false

    std::string curdir = getCurrentDirectory();

    if ( setCurrentDirectory( dir ) ) // if can go to queried dir
    {
        setCurrentDirectory( curdir ); // go back to previous
        return true;
    }

    return false; // otherwise assume dir does not exist
}

bool createDirectory(const std::string& newdir)
{
#ifdef Linux
    // set directory permission, read, write, execute
    // for User, Group and Other
    mode_t perm =   S_IRWXU | S_IRWXG | S_IRWXO;

    if( mkdir( newdir.c_str(), perm ) != 0 )
        return false;
    else
        return true; //success
#endif
}
