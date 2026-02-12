for molfile in $1/*.sdf; do
  sdfile=$(dirname $molfile)/sdf_combined.sdf # remove the .mol extension and add a .sdf
  cat ${molfile} >>${sdfile}
  echo >>${sdfile} # add an extra blank line
  echo '$$$$' >>${sdfile}
done