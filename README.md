## Repository for stellar evolutionary tracks used in TRILEGAL

### Contents

- isotrack/
- isotrack_agb/
- isotrack_parsec/
- isotrack_parcol/

## DOWNLOADING

### IMPORTANT: 

Before proceeding with the next steps DELETE you local copies of: 

- isotrack/
- isortack_agb/
- isotrack_parsec/
- isotrack_parcol/

Your mydirectorytree/trilegal_1.5/.. should now contains:

- trilegal_1.5/

and  if you have already a copy of bolometric corrections tables, also:
 
- photom



The following steps are required to have all the folders in the right place:

1. Move to the trilegal_1.5/.. directory:

   ```
   cd mydirectorytree/trilegal_1.5/..
   ```
   
2. Create an empty Git repository:

   ```
   git init
   ```
   
3. Set the remote "origin":

   ```
   git remote add origin git@gitlab.com:cycyustc/isotracks.git
   ```
   
4. Incorporate changes from the remote repository "isotracks" into the master branch:

   ```
   git pull origin master
   ```

5. To find new commits (i.e. updates in the tracks):

   ```
   git pull origin master
   ```

## VERIFY THAT:

 - mydirectorytree/trilegal_1.5/..  contains the directories

   trilegal_1.5 
   
   isotrack
   
   isotrack_agb
   
   isotrack_parsec

   isotrack_parcol
   
   
   
 - mydirectorytree/trilegal_1.5/ contains the following symbolic link:

   isotrack -> ../isotrack/
   
 - mydirectorytree/trilegal_1.5/isotrack contains the following symbolic links:

   isotrack_agb -> ..isotrack_agb/
   
   parsec -> ../isotrack_parsec/

   parcol -> ../isotrack_parcol/


   
### If you have already a copy of trilegal_1.5 and its input files:

- mydirectorytree/trilegal_1.5/.. should contain the directories

   trilegal_1.5 
   
   isotrack
   
   isotrack_agb
   
   isotrack_parsec

   isotrack_parcol
   
   photom/tab_dust 
   
   photom/bc_dust 
   
   photom/tab_mag_odfnew
   
   photom/bc_odfnew
   
   photom/tab_mag_odfnewbern
   
   photom/bc_odfnewbern
   
   
- mydirectorytree/trilegal_1.5/ instead should contain the following symbolic links:

   isotrack -> ../isotrack/
   
   tab_dust -> ../photom/tab_dust/
   
   bc_dust -> ../photom/bc_dust/
   
   tab_mag_odfnew -> ../photom/tab_mag_odfnew/
   
   bc_odfnew -> ../photom/bc_odfnew/
   
   tab_mag_odfnewbern -> ../photom/tab_mag_odfnewbern/
   
   bc_odfnewbern -> ../photom/bc_odfnewbern/
   
