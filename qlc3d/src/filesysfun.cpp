#include <filesysfun.h>
#include <util/exception.h>
#include <util/logging.h>
namespace FilesysFun
{

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

#ifdef Windows

    bool success = SetCurrentDirectoryA( destdir.c_str() );

    return success;
#endif

}

std::string getCurrentDirectory()
{
    std::string cwd;

#ifdef Linux
    char* buff = get_current_dir_name(); // allocates space to buffer
    cwd = buff;
    free( buff );   // must free buffer
#endif

#ifdef Windows
    char buff[FILENAME_MAX];
    _getcwd( buff, sizeof(buff) );
    cwd = buff;
#endif

    return cwd;

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

#ifdef Windows
    return CreateDirectoryA( newdir.c_str(), NULL );
#endif
}

}//end namespace FilesysFun
