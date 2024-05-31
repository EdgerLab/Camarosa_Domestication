Had to do some hotfixes to panEDTA script, but inconsistent behavior, and would rather not fork it to make these changes.
In future, if having troubles with panEDTA.sh, modify the following:

diff --git a/panEDTA.sh b/panEDTA.sh
index e225dcc..6cde9cb 100755
--- a/panEDTA.sh
+++ b/panEDTA.sh
@@ -113,23 +113,13 @@ echo ""
 
 ## Step 1, initial EDTA annotation, consider to add --sensitive 1, consider to submit each EDTA job to different nodes.
 # make softlink to global cds
-
-# NOTE, SCOTT replaced this
-# if [ $cds != '' ]; then
-#      cds_file=$cds
-#      cds=`basename $cds_file 2>/dev/null`
-#      cds_file=`realpath $cds_file`
-#      ln -s $cds_file $cds 2>/dev/null
-# fi
-
-if [[ "$cds" != '' ]]; then
-    cds_file="$cds"
-    cds=$(basename "$cds_file" 2>/dev/null)
-    cds_file=$(realpath "$cds_file")
-    ln -s "$cds_file" "$cds" 2>/dev/null
+if [ $cds != '' ]; then
+       cds_file=$cds
+       cds=`basename $cds_file 2>/dev/null`
+       cds_file=`realpath $cds_file`
+       ln -s $cds_file $cds 2>/dev/null
 fi

