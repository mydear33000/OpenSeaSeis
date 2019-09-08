/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include "csIOSelection.h"
#include "csFlexHeader.h"
#include "csFlexNumber.h"
#include "csIReader.h"
#include "csException.h"
#include "csSelection.h"
#include "csSortManager.h"
#include "csVector.h"

using namespace std;
using namespace cseis_geolib;

csIOSelection::csIOSelection( std::string const& headerName, int sortOrder, int sortMethod ) {
  myHdrName   = headerName;
  mySortOrder = sortOrder;
  myNumSelectedTraces = 0;
  mySortManager = NULL;
  mySelectedTraceIndexList = new cseis_geolib::csVector<int>();

  if( mySortOrder != csIOSelection::SORT_NONE ) {
    mySortManager = new csSortManager( 1, sortMethod );
  }
  myCurrentSelectedIndex = -1;
  mySelectedValueList = new cseis_geolib::csVector<csFlexNumber*>();
  mySTEPSelection = NULL;
}
csIOSelection::~csIOSelection() {
  if( mySTEPSelection != NULL ) {
    delete mySTEPSelection;
    mySTEPSelection = NULL;
  }
  if( mySelectedValueList != NULL ) {
    for( int i = 0; i < mySelectedValueList->size(); i++ ) {
      delete mySelectedValueList->at(i);
    }
    delete mySelectedValueList;
  }
  if( mySortManager != NULL ) {
    delete mySortManager;
    mySortManager = NULL;
  }
  if( mySelectedTraceIndexList != NULL ) {
    delete mySelectedTraceIndexList;
    mySelectedTraceIndexList = NULL;
  }
}
bool csIOSelection::initialize( cseis_geolib::csIReader* reader, std::string const& hdrValueSelectionText ) {
  step1( reader, hdrValueSelectionText );
  int numTracesToRead = reader->numTraces();
  if( !step2( reader, numTracesToRead ) ) {
    return false;
  }
  return step3( reader );
}


void csIOSelection::step1( cseis_geolib::csIReader* reader, std::string const& hdrValueSelectionText ) {
  cseis_geolib::type_t hdrType;
  reader->setHeaderToPeek( myHdrName, hdrType );
  mySTEPCurrentTraceIndex = reader->getCurrentTraceIndex();  // Save current file trace position
  reader->moveToTrace( 0 );
  mySTEPSelection = new csSelection( 1, &hdrType );
  mySTEPSelection->add( hdrValueSelectionText );
  mySTEPTraceIndex = 0;
}
bool csIOSelection::step2( cseis_geolib::csIReader* reader, int& numTracesToRead ) {
  int lastTrace = std::min( reader->numTraces(), mySTEPTraceIndex+numTracesToRead );
  numTracesToRead = lastTrace - mySTEPTraceIndex;
  for( int itrc = mySTEPTraceIndex; itrc < lastTrace; itrc++ ) {
    csFlexHeader flexHdr;
    if( reader->peekHeaderValue( &flexHdr, itrc ) ) {
      csFlexNumber value(&flexHdr);
      if( mySTEPSelection->contains( &value ) ) {
        mySelectedTraceIndexList->insertEnd(itrc);
        // if( mySortOrder != csIOSelection::SORT_NONE ) 
        mySelectedValueList->insertEnd( new csFlexNumber(&value,mySortOrder == SORT_DECREASING) );
        // Trick to flip polarity of value for sorting in decreasing direction
      }
    }
    else {
      throw( csException("csIOSelection::step2: Error occurred when scanning header value for trace #%d", itrc+1) );
    }
  }
  mySTEPTraceIndex = lastTrace;
  if( mySTEPTraceIndex < reader->numTraces() ) return false;
  return true;
}
bool csIOSelection::step3( cseis_geolib::csIReader* reader ) {
  delete mySTEPSelection;
  mySTEPSelection = NULL;
  // Move reader back to beginning
  if( mySTEPCurrentTraceIndex < 0 ) mySTEPCurrentTraceIndex = 0;
  reader->moveToTrace( mySTEPCurrentTraceIndex );

  // No traces found:
  myNumSelectedTraces = mySelectedTraceIndexList->size();
  if( myNumSelectedTraces == 0 ) {
    return false;
  }

  if( mySortOrder != csIOSelection::SORT_NONE ) {
    mySortManager->resetValues( myNumSelectedTraces );
    for( int is = 0; is < myNumSelectedTraces; is++ ) {
      csFlexNumber* flexNum = mySelectedValueList->at( is );
      //      fprintf(stdout,"Value #%d   %d    %d, trace %d\n", is, flexNum->intValue(), myNumSelectedTraces, mySelectedTraceIndexList->at(is) );
      mySortManager->setValue( is, 0, *flexNum, mySelectedTraceIndexList->at(is) );
    }
    mySelectedTraceIndexList->clear();
    mySortManager->sort();
    if( mySortOrder == SORT_DECREASING ) {
      // Revert polarity flip in sorted values:
      mySortManager->flipPolarity();
    }
  }
  //  else {
  //  for( int is = 0; is < myNumSelectedTraces; is++ ) {
  //    fprintf(stdout,"Value #%d    %d, trace %d\n", is, myNumSelectedTraces, mySelectedTraceIndexList->at(is) );
  //  }
  //}
  myCurrentSelectedIndex = -1;

  return true;
}



int csIOSelection::getCurrentTraceIndex() {
  if( mySortOrder == SORT_NONE ) {
    return mySelectedTraceIndexList->at(myCurrentSelectedIndex);
  }
  else {
    return mySortManager->sortedIndex(myCurrentSelectedIndex);
  }
}
int csIOSelection::getNextTraceIndex() {
  myCurrentSelectedIndex += 1;
  if( myCurrentSelectedIndex >= myNumSelectedTraces ) {
    //throw( csException("csIOSelection::getTraceIndex(): Incorrect call to function. All traces have already been returned.") );
    return -1;
  }
  return getCurrentTraceIndex();
}

int csIOSelection::getNumSelectedTraces() const {
  if( mySortOrder == SORT_NONE ) {
    return mySelectedValueList->size();
  }
  else {
    return mySortManager->numValues();
  }
}

int csIOSelection::getSelectedIndex( int traceIndex ) const {
  if( mySortOrder != SORT_NONE ) {
    return mySortManager->sortedIndex( traceIndex );
  }
  else {
    return traceIndex;
  }
}

cseis_geolib::csFlexNumber const* csIOSelection::getSelectedValue( int traceIndex ) const {
  if( traceIndex >= getNumSelectedTraces() ) {
    throw( csException("csIOSelection::getSelectedValue: Incorrect trace index provided: %d >= %d: This is probably a bug in the calling function",
                       traceIndex, getNumSelectedTraces()) );
  }
  if( mySortOrder == SORT_NONE ) {
    return mySelectedValueList->at( traceIndex );
  }
  else {
    return mySortManager->sortedValue( traceIndex );
  }
}
