.onLoad <- function(lib, pkg)
{
   library.dynam('netresponse', pkg, lib)
   packageStartupMessage('\nnetresponse Copyright (C) 2008-2012 Leo Lahti.\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under GNU GPL >=2, see the licensing terms for details.\n')
}
