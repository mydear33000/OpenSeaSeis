/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_USER_PARAM_H
#define CS_USER_PARAM_H

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "tinyxml2.h"

using namespace tinyxml2;

namespace cseis_system {

/**
* CSeis model list: parsing xml model list file
*
* @author Yuanpeng zhang
* @date   20190909
*/
struct moduleConfig{
    std::string name;
    std::string type;
    int ninport;
    int noutport;
};
class csModuleLists {
public:
  csModuleLists();
  csModuleLists( std::string const& xmlconfigfnm );
  ~csModuleLists();
  void set( std::string const& xmlconfigfnm );

  int const nModules() const;
  std::vector<boost::shared_ptr<moduleConfig>> modulesNMList() const;
  boost::shared_ptr<moduleConfig> getModuleConfig(std::string modelnm);

  std::string const libPath() const;

private:
  std::string modelFullNM(std::string modulenm);

  XMLDocument* xmldoc_;
  std::string libpath_;
  int nmodules_;

  std::vector<boost::shared_ptr<moduleConfig>> moduleconfig_;
};

}

#endif


