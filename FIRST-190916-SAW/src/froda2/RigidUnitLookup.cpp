#include "RigidUnitLookup.h"
#include <algorithm>

void RigidUnitLookup::setLookup_RUtoPlist( 
    const std::vector<IDlist>& RUtoPlist_, size_t nP )
{
  RUtoPlist = RUtoPlist_;
  size_t nRU = RUtoPlist.size();
  
  PtoRUlist.resize( nP );
  PtoRUPlist.resize( nP );
  RUtoRUPlist.resize( nRU );

  //iterate over each rigid unit ru
  int rup = 0;
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    IDlist *plist = &RUtoPlist[ru];
    
    //sort the plist
    sort( plist->begin(), plist->end() );
    
    //iterate over each point p in the plist
    size_t n = plist->size();
    for ( size_t i = 0; i < n; i++ ) {
      int p = (*plist)[i];
      //here we have a point p,
      //a corresponding rigid unit point rup,
      //and its corresponding rigid unit ru.
      //Insert this information in the various lookup tables.
      //Note that the for loop goes over ru in order,
      //and the rup index is incremented in order,
      //so rulists and ruplists will already be sorted as
      //they are built
      PtoRUlist[p].push_back( ru );
      PtoRUPlist[p].push_back( rup );
      RUtoRUPlist[ru].push_back( rup );
      RUPtoRU.push_back( ru );
      RUPtoP.push_back( p );
      rup++;
    }
  }
  
  //here are vector "idioms" for trimming a vector.
  //these lines have the effect of creating a temporary
  //copy of the vector, then copying from the temporary
  //back to the original.  The result is a trimmed vector
  std::vector<IDlist>(PtoRUlist).swap(PtoRUlist);
  std::vector<IDlist>(PtoRUPlist).swap(PtoRUPlist);
  std::vector<int>(RUPtoRU).swap(RUPtoRU);
  std::vector<IDlist>(RUtoRUPlist).swap(RUtoRUPlist);
  std::vector<int>(RUPtoP).swap(RUPtoP);   
}
