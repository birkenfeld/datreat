This is a usefull configuration for xmgrace

to use it copy or link this folder to .grace in your home folder

The most important file is .grace/gracerc as the configguration file
with these lines

DEFINE IFILTER "egrep '^\s*[-+01234456789]|^\s*$' %s" PATTERN "*.dat"

#to define an infutfilter for *.dat files used in Data->Import->ascii
#with  inputfilter "egrep '^\s*[-+01234456789]|^\s*$' %s"  for Pattern "*.dat"

# in templates/defaults.agr
you can change this default file as your standard if wanted


please take a look inside for other stuff