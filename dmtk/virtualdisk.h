#ifndef __DMTK_VDISK_H__
#define __DMTK_VDISK_H__

#include <iostream>
#include <iosfwd>
#include <list>
#include <sstream>
#include <string>
#include <map>
#include <iomanip>


class VirtualDisk:public std::map<std::string,std::string>
{
  public:
    typedef std::map<std::string,std::string> VirtualDiskType;
    typedef std::map<std::string,std::string>::iterator iterator;
    typedef std::map<std::string,std::string>::const_iterator const_iterator;
    VirtualDisk() {}
    
    const std::string& iopen(const char *name); 
    std::string& oopen(const char *name); 

    void close(const std::string &s);

  private:
    iterator pos;
};

std::string &
VirtualDisk::oopen(const char *name)
{
   VirtualDiskType::iterator iter;
   iter = this->find(std::string(name)); 
   if(iter != this->end()) {
      std::string &ss = iter->second;
      ss.clear();
      pos = iter;
      return ss;
   }
   std::string key(name);
   std::string ss;
   std::pair<std::string,std::string> my_pair(key,ss);
   iter = this->insert(this->begin(),my_pair);
   pos = iter;
   return iter->second;
}

const std::string&
VirtualDisk::iopen(const char *name)
{
   VirtualDiskType::iterator iter;
   iter = this->find(std::string(name)); 
   if(iter == this->end()) cout << "ERROR: file " << name << "not found\n";
   pos = iter;
//   if(iter == this->end()) return std::stringstream();
   return iter->second;
}

void
VirtualDisk::close(const std::string &s)
{
  pos->second = s;
}

#endif // __DMTK_VDISK_H__
