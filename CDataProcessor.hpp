#ifndef _DATA_PROCESSOR_
#define _DATA_PROCESSOR_


//Package CImg
#include <CImg.h>

using namespace cimg_library;

#include "CDataBuffer.hpp"

template<typename Tdata, typename Taccess=unsigned char>
class CDataProcessor : public CDataBuffer<Tdata, Taccess>
{
public:
  CImg<Tdata> image;
  //! result access
  CAccessOMPLock laccessR;
  CDataAccess::ACCESS_STATUS_OR_STATE wait_statusR;
  CDataAccess::ACCESS_STATUS_OR_STATE  set_statusR;

  bool do_check;unsigned int check_error;

  CDataProcessor(std::vector<omp_lock_t*> &lock
  , CDataAccess::ACCESS_STATUS_OR_STATE wait_status=CDataAccess::STATUS_FILLED
  , CDataAccess::ACCESS_STATUS_OR_STATE  set_status=CDataAccess::STATUS_PROCESSED
  , CDataAccess::ACCESS_STATUS_OR_STATE wait_statusR=CDataAccess::STATUS_FREE
  , CDataAccess::ACCESS_STATUS_OR_STATE  set_statusR=CDataAccess::STATUS_FILLED
  , bool do_check=false
  )
  : CDataBuffer<Tdata, Taccess>(lock,wait_status,set_status)
  , laccessR(lock[2])
  , wait_statusR(wait_statusR), set_statusR(set_statusR)
  , do_check(do_check)
  {
//    this->debug=true;
    this->class_name="CDataProcessor";

    //this->check_locks(lock);
    if(lock.size()<3)
    {
      printf("code error: locks should have at least 3 locks for %s class.",this->class_name.c_str());
      exit(99);
    }
    check_error=0;
  }//constructor

  //! compution kernel for an iteration
  virtual void kernel(CImg<Tdata> &in,CImg<Tdata> &out)
  {
    std::cout<< __FILE__<<"/"<<__func__<<"(images: in="<<in.width()<<", out="<<out.size()<<") copy kernel, other kernels should be implemented in inherited class."<<std::endl<<std::flush;
    out=in;
  };//kernel

  //! one iteration for any loop
  virtual void iteration(CImg<Taccess> &access,CImgList<Tdata> &images, CImg<Taccess> &accessR,CImgList<Tdata> &results, int n, int i)
  {
    if(this->debug)
    {
      this->lprint.print("",false);
      printf("4 B%02d #%04d: ",n,i);fflush(stdout);
      access.print("access",false);fflush(stderr);
      this->lprint.unset_lock();
    }

    //! 1. compute from buffer
    //wait lock
    unsigned int c=0;
    this->laccess.wait_for_status(access[n],this->wait_status,this->STATE_PROCESSING, c);//filled, processing
    //compution in local
    kernel(images[n],image);

        //check
        if(do_check)
        {
          if(images[n] ==i) NULL; else {++check_error;std::cout<<"compution error: bad generate class for this test."<<std::endl<<std::flush;}
        }

    //unlock
    this->laccess.set_status(access[n],this->STATE_PROCESSING,this->set_status, this->class_name[5],i,n,c);//processing, processed

    //! 2. copy to buffer
    if(this->debug)
    {
      this->lprint.print("",false);
      printf("4 C%02d #%04d: ",n,i);fflush(stdout);
      accessR.print("accessR",false);fflush(stderr);
      this->lprint.unset_lock();
    }

        //check
        if(do_check)
        {
          if(image==i) NULL; else {++check_error;std::cout<<"compution error: bad check (i.e. test failed) on iteration #"<<i<<" (value="<<results[0](0)<<"."<<std::endl<<std::flush;}
        }

    //wait lock
    c=0;
    this->laccessR.wait_for_status(accessR[n],this->wait_statusR,this->STATE_PROCESSING, c);//filled, processing
    //copy local to buffer
    results[n]=image;
    //unlock
    this->laccessR.set_status(accessR[n],this->STATE_PROCESSING,this->set_statusR, this->class_name[5],i,n,c);//processing, processed

  }//iteration

  virtual void show_checking()
  {
    if(do_check)
    {
    if(check_error>0) std::cout<<"test: fail ("<<check_error<<" errors)."<<std::endl;
    else std::cout<<"test: pass."<<std::endl;
    }//do_check
  }//show_checking

};//CDataProcessor

#endif //_DATA_PROCESSOR_

