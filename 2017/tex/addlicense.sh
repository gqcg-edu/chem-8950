for f in *.tex; do 
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" >> tmpfile
  echo "% This work is licensed under the Creative Commons Attribution 4.0 International %" >> tmpfile
  echo "% License. To view a copy of this license, visit                                 %" >> tmpfile
  echo "% http://creativecommons.org/licenses/by/4.0/.                                   %" >> tmpfile
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" >> tmpfile
  cat $f >> tmpfile
  mv tmpfile $f
done
