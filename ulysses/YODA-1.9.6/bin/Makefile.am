dist_bin_SCRIPTS = yoda-config

if ENABLE_PYEXT

## YODA file listing, diffing and modifying
dist_bin_SCRIPTS += yodals yodadiff yodamerge yodastack yodascale

## Histogramming from an ASCII stream of unbinned data
dist_bin_SCRIPTS += yodahist

## Plotting
dist_bin_SCRIPTS += yodacmp yodaplot

## Data format conversion
dist_bin_SCRIPTS += yodacnv yoda2yoda yoda2flat
## Deprecated scripts
dist_bin_SCRIPTS += yoda2aida aida2yoda
## Uber-deprecated scripts
dist_bin_SCRIPTS += aida2flat flat2yoda
if ENABLE_ROOT
dist_bin_SCRIPTS += yoda2root root2yoda
endif

endif

## bash completion
bashcompdir = $(sysconfdir)/bash_completion.d
dist_bashcomp_DATA = yoda-completion
