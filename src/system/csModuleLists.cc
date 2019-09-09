/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csModuleLists.h"
#include "cseis_defines.h"
#include "csException.h"
#include "geolib_string_utils.h"
#include <boost/filesystem.hpp>
#include "glog/logging.h"

using namespace cseis_system;
namespace fs = boost::filesystem;

csModuleLists::csModuleLists()
    : xmldoc_(nullptr)
    , nmodules_(0)
{}
csModuleLists::csModuleLists( std::string const& xmlconfigfnm )
    : xmldoc_(nullptr)
    , nmodules_(0)
{
    set(xmlconfigfnm);
}
csModuleLists::~csModuleLists(){
    moduleconfig_.clear();
    if(xmldoc_) delete xmldoc_;
}
void csModuleLists::set( std::string const& xmlconfigfnm ){
    if(xmldoc_) delete xmldoc_;
    nmodules_ = 0;

    moduleconfig_.clear();

    xmldoc_ = new XMLDocument();
    XMLError error = xmldoc_->LoadFile(xmlconfigfnm.c_str());
    if (error != XML_NO_ERROR)
        LOG(ERROR) << "Error loading config file " << xmlconfigfnm;
    XMLElement* configTag = xmldoc_->FirstChildElement("config");
    libpath_.assign(configTag->FirstChildElement("libdir")->GetText());
    if(libpath_.at(libpath_.length() - 1) != FORWARD_SLASH)libpath_.append("/");

    // check lib exist!
    if(!fs::is_directory(libpath_))
        LOG(ERROR) << "Can not find lib path " << libpath_;

    XMLElement* moduleTag = configTag->FirstChildElement("module");
    while(moduleTag){
        boost::shared_ptr<moduleConfig> m(new moduleConfig);
        m->name.assign(moduleTag->Attribute("name"));
        m->type.assign(moduleTag->Attribute("type"));
        m->ninport = moduleTag->IntAttribute("inport");
        m->noutport = moduleTag->IntAttribute("outport");
        CHECK_GE(m->ninport, 0) << "Module " << m->name << " inport " << m->ninport << " < 0?";
        CHECK_LE(m->ninport, 2) << "Module " << m->name << " inport " << m->ninport << " > 2?";
        CHECK_GE(m->noutport, 0) << "Module " << m->name << " outport " << m->noutport << " < 0?";
        CHECK_LE(m->noutport, 2) << "Module " << m->name << " outport " << m->noutport << " > 2?";

        if(!fs::exists(modelFullNM(m->name)))
            LOG(ERROR) << "Can't find module " << m->name;

        moduleconfig_.push_back(m);
        nmodules_++;

        moduleTag = moduleTag->NextSiblingElement("module");
    }
}
std::string csModuleLists::modelFullNM(std::string modulenm){
    std::string libnm = libpath_ + MODULEPREFIX + modulenm + ".so";

    return libnm;
}
int const csModuleLists::nModules() const{
    return nmodules_;
}
std::vector<boost::shared_ptr<moduleConfig>> csModuleLists::modulesNMList() const{
    return moduleconfig_;
}
boost::shared_ptr<moduleConfig> csModuleLists::getModuleConfig(std::string modelnm){    
    for(int i = 0; i != nmodules_; i++){
        if(!moduleconfig_.at(i)->name.compare(modelnm))
            return moduleconfig_.at(i);
    }
    return nullptr;
}
std::string const csModuleLists::libPath() const{
    return libpath_;
}


