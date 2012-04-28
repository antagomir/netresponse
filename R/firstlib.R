# Copyright (C) 2008-2012 Olli-Pekka Huovilainen and Leo Lahti 
# Contact: Leo Lahti <leo.lahti@iki.fi>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

.onLoad <- function(lib, pkg)
{
   library.dynam('netresponse', pkg, lib)
   packageStartupMessage('\nnetresponse Copyright (C) 2008-2012 Leo Lahti.\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under GNU GPL >=2, see the licensing terms for details.\n')
}
