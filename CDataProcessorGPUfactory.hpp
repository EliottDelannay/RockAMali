#ifndef _DATA_PROCESSOR_GPU_FACTORY_
#define _DATA_PROCESSOR_GPU_FACTORY_

#include "CDataProcessorGPU.hpp"

//! factory for GPU processing
/**
 *
**/
template<typename Tdata, typename Taccess=unsigned char>
class CDataProcessorGPUfactory
{
public:
  static CDataProcessorGPU<Tdata, Taccess> *NewCDataProcessorGPU(const std::string &name//="list types"
  , std::vector<std::string> &factory_types
  , std::vector<omp_lock_t*> &lock
  , compute::device device, int VECTOR_SIZE
  , CDataAccess::ACCESS_STATUS_OR_STATE wait_status=CDataAccess::STATUS_FILLED
  , CDataAccess::ACCESS_STATUS_OR_STATE  set_status=CDataAccess::STATUS_PROCESSED
  , CDataAccess::ACCESS_STATUS_OR_STATE wait_statusR=CDataAccess::STATUS_FREE
  , CDataAccess::ACCESS_STATUS_OR_STATE  set_statusR=CDataAccess::STATUS_FILLED
  , bool do_check=false
  )
  {
    //reset
    factory_types.clear();
    //if
    factory_types.push_back      ("copy")         ;if(name == factory_types.back())
      return new CDataProcessorGPU<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);
    factory_types.push_back      ("program")             ;if(name == factory_types.back())
      return new CDataProcessorGPU_opencl<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);
    factory_types.push_back      ("lambda")              ;if(name == factory_types.back())
      return new CDataProcessorGPU_lambda<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);
    factory_types.push_back      ("closure")              ;if(name == factory_types.back())
      return new CDataProcessorGPU_closure<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);
    factory_types.push_back      ("function")              ;if(name == factory_types.back())
      return new CDataProcessorGPU_function<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);
    factory_types.push_back      ("function_lambda")              ;if(name == factory_types.back())
      return new CDataProcessorGPU_function_lambda<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);
    factory_types.push_back      ("function_macro")              ;if(name == factory_types.back())
      return new CDataProcessorGPU_function_macro<Tdata, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check);

//    if(name == "error")
//      return new CDataProcessorGPU_error;
    //listing known types in factory
    if(name=="list types")
      return NULL;
    //unknown
    std::cerr<<"Module name is unknown, i.e. \""<<name<<"\"."<<std::endl;
    return NULL;
  }//NewCDataProcessorGPU

  //! get type list in factory
  static void get_factory_types(std::vector<std::string> &factory_types)
  {
    //dummy vars
    std::vector<omp_lock_t*> locks;
    ///dummy OpenMP locks
    omp_lock_t print_lock;omp_init_lock(&print_lock);
    locks.push_back(&print_lock);
    ///dummy OpenCL device
    compute::device device=compute::device();//=compute::system::default_device();
    CDataProcessorGPUfactory<Tdata, Taccess>::NewCDataProcessorGPU("list types",factory_types /*dummies then default*/,locks,device,0);
  }//get_factory_types
  //! show type list in factory
  static void show_factory_types(std::vector<std::string> &factory_types)
  {
    if(factory_types.empty()) get_factory_types(factory_types);
    std::cout<<"information: -factory- GPU processing types are: ";
    unsigned int i;
    for(i=0;i<factory_types.size()-1;++i)
      std::cout<<factory_types[i]<<", ";
    std::cout  <<factory_types[i]<<".";
  }//show_factory_types
};//CDataProcessorGPUfactory

#endif //_DATA_PROCESSOR_GPU_FACTORY_

