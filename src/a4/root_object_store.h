#ifndef _A4_STORE_ROOT_OBJECT_STORE_H_
#define _A4_STORE_ROOT_OBJECT_STORE_H_

#ifdef HAVE_CERN_ROOT_SYSTEM

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <utility>

#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TAxis.h>


#include <a4/storable.h>
#include <a4/object_store.h>

namespace a4 {
namespace store {

class RootStorable : public Storable {
public:
    virtual TObject* get_tobject() = 0;
};


// TODO: Figure out how different dimensionalities of histograms will work

template<class ROOT_TYPE>
class SpecificRootStorable : public RootStorable {
public:
    typedef SpecificRootStorable<ROOT_TYPE> This;

    SpecificRootStorable() : _initializations_remaining(1) {}


#ifdef A4STORE_STANDALONE
            // TODO: implement alternative storage mechanisms
#else
    virtual shared<const google::protobuf::Message> as_message() { FATAL("Can't serialize root objects to protobuf"); }
    virtual void construct_from(const google::protobuf::Message&) { FATAL("Can't build root obects from protobuf"); }
    virtual void construct_from(shared<google::protobuf::Message> m) { construct_from(*m); }
#endif // ifdef A4STORE_STANDALONE

    virtual Storable& operator+=(const Storable &other) { abort(); }
    virtual Storable&& operator+(const Storable &other) { abort(); }
    virtual Storable& operator*=(const double&) { abort(); }
    virtual Storable&& clone_storable() { abort(); }
    
    ROOT_TYPE root_object;
    
    TObject* get_tobject() { return &root_object; }


    virtual void constructor(const char* title) {
        _initializations_remaining++;
        root_object.SetTitle(title);
    }
    virtual void constructor(const uint32_t bins, const double min, const double max, const char* label="") {
      FATAL("Constructor not implemented");
    }
    virtual void constructor(const std::vector<double>& bins, const char* label="") {
      FATAL("Constructor not implemented");
    }
    
    template <typename... Args> This& operator()(const Args&... args) {
        if (_initializations_remaining != 0) {
            _initializations_remaining--;
            static_cast<This*>(this)->constructor(args...);
        }
        return *static_cast<This*>(this);
    }
    
    int _initializations_remaining;
    bool _initialized() const { return _initializations_remaining == 0; }
    
    virtual void fill(double v, double w=1) {
      FATAL("Fill not implemented");
    }
};

  typedef SpecificRootStorable<TH1D> RTH1;

template<class ROOT_TYPE, int DIMENSIONS>
class StorableRootHist : public SpecificRootStorable<ROOT_TYPE> {
public:

  StorableRootHist() {this->_initializations_remaining = DIMENSIONS;}


  void constructor(const std::vector<double>& bins, const char* label="") {
    TAxis* current_axis(0);
    int n_cells(0);

    unsigned int nbins(bins.size()-1);
    assert(nbins > 0);

    // Set the correct axis/cells based on the number of dimensions and itrerations remaining
    if (this->_initializations_remaining == DIMENSIONS - 3) {    
      current_axis = this->root_object.GetZaxis();
      n_cells = (nbins + 2) * (this->root_object.GetNbinsY() + 2) * (this->root_object.GetNbinsX() + 2); 
    } else if (this->_initializations_remaining == DIMENSIONS - 2) {
      current_axis = this->root_object.GetYaxis();
      n_cells = (nbins + 2) * (this->root_object.GetNbinsX() + 2); 
    } else if (this->_initializations_remaining == DIMENSIONS - 1) {
      current_axis = this->root_object.GetXaxis();
      n_cells = nbins + 2;
    }

    // Init the current axis 
    current_axis->SetRange(0,0);
    current_axis->Set(nbins, &bins[0]);
    current_axis->SetTitle(label);

    // Update the histogram
    this->root_object.SetBinsLength(n_cells);

    if (this->root_object.GetSumw2N()) {
      this->root_object.GetSumw2()->Set(n_cells);
    }

    this->root_object.Reset();
  }

  void constructor(const uint32_t bins, const double min, const double max, const char* label="") {
    TAxis* current_axis(0);
    int n_cells(0);

    // Set the correct axis/cells based on the number of dimensions and itrerations remaining
    if (this->_initializations_remaining == DIMENSIONS - 3) {    
      current_axis = this->root_object.GetZaxis();
      n_cells = (bins + 2) * (this->root_object.GetNbinsY() + 2) * (this->root_object.GetNbinsX() + 2); 
    } else if (this->_initializations_remaining == DIMENSIONS - 2) {
      current_axis = this->root_object.GetYaxis();
      n_cells = (bins + 2) * (this->root_object.GetNbinsX() + 2); 
    } else if (this->_initializations_remaining == DIMENSIONS - 1) {
      current_axis = this->root_object.GetXaxis();
      n_cells = bins + 2;
    }

    // Init the current axis 
    current_axis->SetRange(0,0);
    current_axis->Set(bins, min, max);
    current_axis->SetTitle(label);

    // Update the histogram
    this->root_object.SetBinsLength(n_cells);

    if (this->root_object.GetSumw2N()) {
      this->root_object.GetSumw2()->Set(n_cells);
    }

    this->root_object.Reset();
  }

  void fill(double x, double w=1) {
    //FATAL("Fill not implemented");
     assert(this->_initialized());
     this->root_object.Fill(x, w);
  }
  
};

// class RTH1D : public StorableRootHist<ROOT_TYPE, 1> {
// public:
//   
//   RTH1D() {this->StorableRootHist();}
// 
//   void fill(double x, double w=1) {
//     assert(this->_initialized());
//     this->root_object.Fill(x, w);
//   }
// 
// };


  typedef StorableRootHist<TH1D, 1> RTH1D;
  typedef StorableRootHist<TH2D, 2> RTH2D;
  typedef StorableRootHist<TH3D, 2> RTH3D;

class RootObjectStore : public ObjectBackStore {
public:

    static TDirectory* mkdirs(TDirectory* start, const std::string& file) {
        if (!start->GetDirectory(file.c_str()))
            start->mkdir(file.c_str());
        TDirectory* d = start->GetDirectory(file.c_str());
        d->cd();
        return d;
    }
    
  static std::pair<std::string, std::string> path_and_name(const std::string& full) {
    std::string::size_type idelim = full.find_last_of("/");

    std::string path, name;

    if (idelim == std::string::npos) {
      path = "";
      name = full;     
    } else {
      path = full.substr(0, idelim);      
      name = full.substr(idelim+1);      
    }

    return std::make_pair(path, name);
  }


  void write(std::string filename = "") {
        TDirectory* destination = gDirectory;
        
	if (!filename.empty()) {
	  destination = new TFile(filename.c_str(), "RECREATE");
	}

        assert(destination->IsWritable());
        
        for (std::map<std::string, shared<Storable> >::iterator i = _store->begin(); 
            i != _store->end(); i++) {

            RootStorable* object = dynamic_cast<RootStorable*>(i->second.get());

            if (!object) {
                std::cout << "Non-RootStorable object: " << i->first << std::endl;
                continue;
            }
            
	    std::pair<std::string, std::string> path = path_and_name(i->first);
            TDirectory* d = mkdirs(destination, path.first);
            d->WriteTObject(object->get_tobject(), path.second.c_str());
        }
        
    }
};
    

}
}
    
#endif // HAVE_CERN_ROOT_SYSTEM

#endif // _A4_STORE_ROOT_OBJECT_STORE_H_

