#ifndef _DATA_PROCESSOR_GPU_OPENCL_
#define _DATA_PROCESSOR_GPU_OPENCL_

#include "CDataProcessorGPU.hpp"

//! complex operation with OpenCL including template types for GPU process
/**
 *  FMA: val * 2.1 + 123.45
 *  \note: Tdata and Tproc could be in the template source
**/
template<typename Tdata=unsigned int,typename Tproc=unsigned int, typename Taccess=unsigned char>
class CDataProcessorGPU_opencl_template : public CDataProcessorGPU_vMcPc_check<Tdata,Tproc, Taccess>
{
public:
  compute::program program;
  compute::kernel  kernel;
  bool kernel_loaded;
//OpenCL function for this class
  //! OpenCL kernel name
  std::string kernel_name;
  //! OpenCL source string (with template)
  std::string source_with_template;
//! OpenCL source (with template)
/**
 * template types must be either \c Tdata or \c Tproc
 * \note this function should redefined in inherited class
**/
void define_opencl_source()
{
  kernel_name="vMcPc";
  source_with_template=BOOST_COMPUTE_STRINGIZE_SOURCE(
  __kernel void vMcPc(__global const Tdata*input, int size, __global Tproc*output)
  {
    const int gid = get_global_id(0);
    output[gid]=input[gid]*2.1+123.45;
  }
  );//source with template
}//define_opencl_source
//! OpenCL program from source (with template)
/**
 * template types: \c Tdata or \c Tproc only will be tranlated, not more !
**/
compute::program make_opencl_program(const compute::context& context)
{
  //init. source
  define_opencl_source();
  //translate template
  std::string source=source_with_template;
  std::vector<std::string> str_old;str_old.push_back(    "Tdata");              str_old.push_back(    "Tproc");
  std::vector<std::string> str_new;str_new.push_back(CImg<Tdata>::pixel_type());str_new.push_back(CImg<Tproc>::pixel_type());
  //! \bug [template in OpenCL] CImg<T>::pixel_type() is almost ok, except for "long int" giving "int64" (and may be more !?)
  for(unsigned int i=0;i<str_old.size();++i)
  {
    //replace all str_old by str_new
    std::string::size_type pos=0;
    while( (pos=source.find(str_old[i], pos)) != std::string::npos )
    {
      source.replace(pos,str_old[i].size(), str_new[i]);
      pos+=str_new[i].size();
    }//replace loop
  }//for loop
std::cout<<"source:"<<std::endl<<"\""<<source<<std::endl<<"\""<<std::endl<<std::flush;
  // create program
  return compute::program::build_with_source(source,context);
}//make_opencl_program

  CDataProcessorGPU_opencl_template(std::vector<omp_lock_t*> &lock
  , compute::device device, int VECTOR_SIZE
  , CDataAccess::ACCESS_STATUS_OR_STATE wait_status=CDataAccess::STATUS_FILLED
  , CDataAccess::ACCESS_STATUS_OR_STATE  set_status=CDataAccess::STATUS_PROCESSED
  , CDataAccess::ACCESS_STATUS_OR_STATE wait_statusR=CDataAccess::STATUS_FREE
  , CDataAccess::ACCESS_STATUS_OR_STATE  set_statusR=CDataAccess::STATUS_FILLED
  , bool do_check=false
  )
  : CDataProcessorGPU_vMcPc_check<Tdata,Tproc, Taccess>(lock,device,VECTOR_SIZE,wait_status,set_status,wait_statusR,set_statusR,do_check)
  {
    this->debug=true;
    this->check_locks(lock);
    //OpenCL framework
    program=make_opencl_program(this->ctx);
    this->class_name="CDataProcessorGPU_openclT_"+kernel_name;
    kernel_loaded=false;
  }//constructor

  //! compution kernel for an iteration (compution=copy, here)
  virtual void kernelGPU(compute::vector<Tdata> &in,compute::vector<Tproc> &out)
  {
    if(!kernel_loaded)
    {//load kernel
      kernel=compute::kernel(program,kernel_name.c_str());
      kernel.set_arg(0,this->device_vector_in.get_buffer());
      kernel.set_arg(1,(int)this->device_vector_in.size());
      kernel.set_arg(2,this->device_vector_out.get_buffer());
      kernel_loaded=true;
    }//load kernel once
    //compute
    using compute::uint_;
    uint_ tpb=16;
    uint_ workSize=this->device_vector_in.size();
    this->queue.enqueue_1d_range_kernel(kernel,0,workSize,tpb);
  };//kernelGPU

};//CDataProcessorGPU_opencl_template

#endif //_DATA_PROCESSOR_GPU_OPENCL_

