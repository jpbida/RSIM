use 5.008002;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'bida::pdb::chain',
    VERSION_FROM      => 'lib/bida/pdb/chain.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/bida/pdb/chain.pm', # retrieve abstract from module
       AUTHOR         => 'John Bida;bida.john@mayo.edu;14489813 <m003206@mayo.edu>') : ()),
);
