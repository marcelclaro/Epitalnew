/*
 * TYPE.cpp
 *
 *  Created on: Jul 31, 2012
 *      Author: marcel
 */


#include <vector>
#include <list>
#include <stdexcept>
#include <iostream>
#include <typeinfo>
#include "Heterostructure.hpp"

#ifndef Heterostructure_CPP_
#define Heterostructure_CPP_

namespace epital{

template <typename TYPE>
Heterostructure<TYPE>::Epilayer::Epilayer(shared_ptr<Material> material, TYPE layerwidth, TYPE doping_): width(layerwidth), doping(doping_) {
	if(material)
		mat=material;
	else
		throw invalid_argument("Invalid 'shared_ptr<Material>': Pointing to null Material");
}

template <typename TYPE>
Heterostructure<TYPE>::Epilayer::~Epilayer() {

}

template <typename TYPE>
void Heterostructure<TYPE>::Epilayer::setTemperature(Temperature temp){
	mat->setTemperature(temp);
}




template <typename TYPE>
shared_ptr<Material> Heterostructure<TYPE>::Epilayer::getMaterial(){
	return mat;
}

template <typename TYPE>
TYPE Heterostructure<TYPE>::Epilayer::getWidth() const{
	return width;
}

template <typename TYPE>
TYPE Heterostructure<TYPE>::Epilayer::getDoping() const{
	return doping;
}

template <typename TYPE>
Heterostructure<TYPE>::Heterostructure(vector<Epilayer>& activeregion_, vector<Epilayer>& contacts_, shared_ptr<Material> substrate_, int roughness): activeregion(activeregion_), contacts(contacts_), substrate(substrate_), bias(0), ladderbias(0), hascontact(false) {
	if(contacts.size()==2){
		hascontact=true;
		TYPE currpos = -(contacts[0].getWidth());
		limitslst.push_back(make_pair(currpos,0.0));
		currpos=0;
		matlst.push_back(contacts[0].getMaterial());
		widthlst.push_back(contacts[0].getWidth());
		dopinglst.push_back(contacts[0].getDoping());

		if(roughness!=0){
			vector<Epilayer> newactiveregion;
			for(auto lst=activeregion.begin();lst!=(activeregion.end()-1);lst++){
				TYPE monolayerwidth=lst->getMaterial()->latticeparam()/2.0;
				newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(lst->getMaterial(),lst->getWidth()-(std::ceil((roughness-1)/2.0)+0.5)*monolayerwidth,lst->getDoping()));
				std::list<shared_ptr<Material> > appendlstbefore,appendlstafter;
				std::list<TYPE> appenddopinglst,appenddopinglstafter;
				appendlstbefore.push_back(mixMaterial(lst->getMaterial(),(lst+1)->getMaterial()));
				appendlstafter.push_back(mixMaterial(lst->getMaterial(),(lst+1)->getMaterial()));
				appenddopinglst.push_back((lst->getDoping()+(lst+1)->getDoping())/2.0);
				appenddopinglstafter.push_back((lst->getDoping()+(lst+1)->getDoping())/2.0);
				for(int a=std::ceil((roughness-1)/2.0);a!=0;a--){
					appendlstbefore.push_back(mixMaterial(lst->getMaterial(),appendlstbefore.back()));
					appendlstafter.push_back(mixMaterial((lst+1)->getMaterial(),appendlstafter.back()));
					appenddopinglst.push_back((lst->getDoping()+appenddopinglst.back())/2.0);
					appenddopinglstafter.push_back(((lst+1)->getDoping()+appenddopinglstafter.back())/2.0);
				}
				if(roughness-1>0){
					for(int b=std::ceil((roughness-1)/2.0);b>0;b--){
						newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(appendlstbefore.back(),monolayerwidth,appenddopinglst.back()));
						appendlstbefore.pop_back();
						appenddopinglst.pop_back();
					}
					for(int c=std::ceil((roughness-1)/2.0);c>=0;c--){
						newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(appendlstafter.front(),monolayerwidth,appenddopinglstafter.front()));
						appendlstafter.pop_front();
						appenddopinglstafter.pop_front();
					}

				}
				else
					newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(appendlstafter.front(),monolayerwidth,appenddopinglstafter.front()));

				(lst+1)->width-=(std::ceil((roughness-1)/2.0)+0.5)*monolayerwidth;
			}
			newactiveregion.push_back(Heterostructure<TYPE>::Epilayer((activeregion.end()-1)->getMaterial(),(activeregion.end()-1)->getWidth(),(activeregion.end()-1)->getDoping()));
			activeregion=newactiveregion;
		}

		for(auto lst : activeregion)
			matlst.push_back(lst.getMaterial());

		for(auto lst : activeregion)
			widthlst.push_back(lst.getWidth());

		for(auto lst : activeregion)
			dopinglst.push_back(lst.getDoping());

		for(auto lst : activeregion){
			limitslst.push_back(make_pair(currpos,currpos+lst.getWidth()));
			currpos+=lst.getWidth();
		}



		matlst.push_back(contacts[1].getMaterial());
		widthlst.push_back(contacts[1].getWidth());
		dopinglst.push_back(contacts[1].getDoping());
		limitslst.push_back(make_pair(currpos,currpos+contacts[1].getWidth()));


		TYPE exx_yy;
		TYPE ezz;
		TYPE normalstrain;
		TYPE qfactorstrain;
		TYPE p;
		TYPE q;
		for(auto lst : matlst){
			exx_yy=(substrate->latticeparam()-lst->latticeparam())/lst->latticeparam();
			ezz=-2*(lst->getC12()/lst->getC11())*exx_yy;
			normalstrain=(2.0*exx_yy+ezz);
			qfactorstrain=(2.0*exx_yy-2.0*ezz);
			p=(lst->getStrain_av()*normalstrain);
			q=-(lst->getStrain_b()*0.5*qfactorstrain);
			dEc.push_back(lst->getStrain_ac()*normalstrain);
			dEv_HH.push_back(-p-q);
			dEv_lH.push_back(-p+q);
		}

		setBIAS(bias);
		setBIASladder(ladderbias,activeregion.size());

	}
	else if(contacts.size()==0){
		TYPE currpos=0;
		if(roughness!=0){
			vector<Epilayer> newactiveregion;
			for(auto lst=activeregion.begin();lst!=(activeregion.end()-1);lst++){
				TYPE monolayerwidth=lst->getMaterial()->latticeparam()/2.0;
				newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(lst->getMaterial(),lst->getWidth()-(std::ceil((roughness-1)/2.0)+0.5)*monolayerwidth,lst->getDoping()));
				std::list<shared_ptr<Material> > appendlstbefore,appendlstafter;
				std::list<TYPE> appenddopinglst,appenddopinglstafter;
				appendlstbefore.push_back(mixMaterial(lst->getMaterial(),(lst+1)->getMaterial()));
				appendlstafter.push_back(mixMaterial(lst->getMaterial(),(lst+1)->getMaterial()));
				appenddopinglst.push_back((lst->getDoping()+(lst+1)->getDoping())/2.0);
				appenddopinglstafter.push_back((lst->getDoping()+(lst+1)->getDoping())/2.0);
				for(int a=std::ceil((roughness-1)/2.0);a!=0;a--){
					appendlstbefore.push_back(mixMaterial(lst->getMaterial(),appendlstbefore.back()));
					appendlstafter.push_back(mixMaterial((lst+1)->getMaterial(),appendlstafter.back()));
					appenddopinglst.push_back((lst->getDoping()+appenddopinglst.back())/2.0);
					appenddopinglstafter.push_back(((lst+1)->getDoping()+appenddopinglstafter.back())/2.0);
				}
				if(roughness-1>0){
					for(int b=std::ceil((roughness-1)/2.0);b>0;b--){
						newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(appendlstbefore.back(),monolayerwidth,appenddopinglst.back()));
						appendlstbefore.pop_back();
						appenddopinglst.pop_back();
					}
					for(int c=std::ceil((roughness-1)/2.0);c>=0;c--){
						newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(appendlstafter.front(),monolayerwidth,appenddopinglstafter.front()));
						appendlstafter.pop_front();
						appenddopinglstafter.pop_front();
					}

				}
				else
					newactiveregion.push_back(Heterostructure<TYPE>::Epilayer(appendlstafter.front(),monolayerwidth,appenddopinglstafter.front()));
				(lst+1)->width-=(std::ceil((roughness-1)/2.0)+0.5)*monolayerwidth;
			}
			newactiveregion.push_back(Heterostructure<TYPE>::Epilayer((activeregion.end()-1)->getMaterial(),(activeregion.end()-1)->getWidth(),(activeregion.end()-1)->getDoping()));
			activeregion=newactiveregion;
		}

		for(auto lst : activeregion)
			matlst.push_back(lst.getMaterial());

		for(auto lst : activeregion)
			widthlst.push_back(lst.getWidth());

		for(auto lst : activeregion)
			dopinglst.push_back(lst.getDoping());

		for(auto lst : activeregion){
			limitslst.push_back(make_pair(currpos,currpos+lst.getWidth()));
			currpos+=lst.getWidth();
		}

		TYPE exx_yy;
		TYPE ezz;
		TYPE normalstrain;
		TYPE qfactorstrain;
		TYPE p;
		TYPE q;
		for(auto lst : matlst){
			exx_yy=(substrate->latticeparam()-lst->latticeparam())/lst->latticeparam();
			ezz=-2*(lst->getC12()/lst->getC11())*exx_yy;
			normalstrain=(2.0*exx_yy+ezz);
			qfactorstrain=(2.0*exx_yy-2.0*ezz);
			p=-(lst->getStrain_av()*normalstrain);
			q=-(lst->getStrain_b()*0.5*qfactorstrain);
			dEc.push_back(lst->getStrain_ac()*normalstrain);
			dEv_HH.push_back(-p-q);
			dEv_lH.push_back(-p+q);
		}

		setBIAS(bias);
		setBIASladder(ladderbias,matlst.size());

	}
	else{
		throw runtime_error("Invalid number of contacts");
	}
}

template <typename TYPE>
Heterostructure<TYPE>::~Heterostructure() {

}

template <typename TYPE>
TYPE Heterostructure<TYPE>::getTotalDoping(){
	long int currlayer=0;
	TYPE sum=0;
	for(auto i : dopinglst){
		sum+=i*widthlst[currlayer];
		++currlayer;
	}
	return sum;
}

template <typename TYPE>
bool Heterostructure<TYPE>::hasContact(){
	return hascontact;
}

template <typename TYPE>
shared_ptr<Material> Heterostructure<TYPE>::getSubstrate(){
	return substrate;
}


template <typename TYPE>
void Heterostructure<TYPE>::duplicate() {
	if(contacts.size()==2){
		matlst.pop_back();
		for(auto lst : activeregion)
			matlst.push_back(lst.getMaterial());
		matlst.push_back(contacts[1].getMaterial());

		widthlst.pop_back();
		for(auto lst : activeregion)
			widthlst.push_back(lst.getWidth());
		widthlst.push_back(contacts[1].getWidth());

		dopinglst.pop_back();
		for(auto lst : activeregion)
			dopinglst.push_back(lst.getDoping());
		dopinglst.push_back(contacts[1].getDoping());

		TYPE currpos = -(contacts[0].getWidth());
		limitslst.pop_back();
		currpos=0;
		for(auto lst : activeregion){
			limitslst.push_back(make_pair(currpos,currpos+lst.getWidth()));
			currpos+=lst.getWidth();
		}
		limitslst.push_back(make_pair(currpos,currpos+contacts[1].getWidth()));

		TYPE exx_yy;
		TYPE ezz;
		TYPE normalstrain;
		TYPE qfactorstrain;
		TYPE p;
		TYPE q;
		for(auto lst : matlst){
			exx_yy=(substrate->latticeparam()-lst->latticeparam())/lst->latticeparam();
			ezz=-2*(lst->getC12()/lst->getC11())*exx_yy;
			normalstrain=(2.0*exx_yy+ezz);
			qfactorstrain=(2.0*exx_yy-2.0*ezz);
			p=(lst->getStrain_av()*normalstrain);
			q=-(lst->getStrain_b()*0.5*qfactorstrain);
			dEc.push_back(lst->getStrain_ac()*normalstrain);
			dEv_HH.push_back(-p-q);
			dEv_lH.push_back(-p+q);
		}

		auto activelayercopy = activeregion;
		for(auto k : activelayercopy){
			activeregion.push_back(k);
		}

		setBIAS(2*bias);
		setBIASladder(2*ladderbias,activeregion.size());

	}
	else if(contacts.size()==0){
		for(auto lst : activeregion)
			matlst.push_back(lst.getMaterial());

		for(auto lst : activeregion)
			widthlst.push_back(lst.getWidth());

		for(auto lst : activeregion)
			dopinglst.push_back(lst.getDoping());

		TYPE currpos = limitslst.back().second;
		for(auto lst : activeregion){
			limitslst.push_back(make_pair(currpos,currpos+lst.getWidth()));
			currpos+=lst.getWidth();
		}

		TYPE exx_yy;
		TYPE ezz;
		TYPE normalstrain;
		TYPE qfactorstrain;
		TYPE p;
		TYPE q;
		for(auto lst : matlst){
			exx_yy=(substrate->latticeparam()-lst->latticeparam())/lst->latticeparam();
			ezz=-2*(lst->getC12()/lst->getC11())*exx_yy;
			normalstrain=(2.0*exx_yy+ezz);
			qfactorstrain=(2.0*exx_yy-2.0*ezz);
			p=-(lst->getStrain_av()*normalstrain);
			q=-(lst->getStrain_b()*0.5*qfactorstrain);
			dEc.push_back(lst->getStrain_ac()*normalstrain);
			dEv_HH.push_back(-p-q);
			dEv_lH.push_back(-p+q);
		}

		auto activelayercopy = activeregion;
		for(auto k : activelayercopy){
			activeregion.push_back(k);
		}

		setBIAS(2*bias);
		setBIASladder(2*ladderbias,activeregion.size());

	}
	else{
		throw runtime_error("Invalid number of contacts");
	}
}


template <typename TYPE>
void Heterostructure<TYPE>::setBIAS(TYPE BIAS){
	bias=BIAS;

	TYPE activewidth = 0;
	TYPE accumulated = 0;
	TYPE electricfield;

	linearelectricfield.clear();
	layerbeginenergy.clear();

	if(hascontact){

		long int index = 0;
		for(TYPE i : widthlst){
			activewidth+=i/matlst[index]->getDielectric_static();
			++index;
		}

		activewidth-= contacts[0].getWidth() / (matlst[index])->getDielectric_static();
		activewidth-= contacts[1].getWidth() / (matlst[index])->getDielectric_static();

		linearelectricfield.push_back(0.0);
		layerbeginenergy.push_back(accumulated);
		for(Epilayer j : activeregion){
			electricfield = -bias/(j.getMaterial()->getDielectric_static()*activewidth);
			linearelectricfield.push_back(electricfield);
			accumulated += Constant::e*electricfield*j.getWidth();
			layerbeginenergy.push_back(accumulated);
		}
		linearelectricfield.push_back(0.0);
		layerbeginenergy.push_back(accumulated);
	}
	else{

		long int index = 0;
		for(TYPE i : widthlst){
			activewidth+=i/matlst[index]->getDielectric_static();
			++index;
		}

		layerbeginenergy.push_back(accumulated);
		for(Epilayer j : activeregion){
			electricfield = -bias/(j.getMaterial()->getDielectric_static()*activewidth);
			linearelectricfield.push_back(electricfield);
			accumulated += Constant::e*electricfield*j.getWidth();
			layerbeginenergy.push_back(accumulated);
		}
		layerbeginenergy.push_back(accumulated);
	}

}


template <typename TYPE>
void Heterostructure<TYPE>::setBIASladder(TYPE BIAS, long int layersperiod){
	ladderbias = BIAS;

	electricfieldladder.clear();

	long int periods = activeregion.size()/layersperiod;

	auto biasperperiod = ladderbias/periods;
	TYPE currentbias = 0;

	if(hascontact){
		electricfieldladder.push_back(0);

		for(auto i : activeregion){
			for(long int k=0; k < layersperiod; k++)
				electricfieldladder.push_back(-Constant::e*currentbias);
			currentbias+=biasperperiod;
		}

		electricfieldladder.push_back(currentbias-biasperperiod);

	}
	else{
			for(auto i : activeregion){
				for(long int k=0; k < layersperiod; k++)
					electricfieldladder.push_back(-Constant::e*currentbias);
				currentbias+=biasperperiod;
			}

	}

}

template <typename TYPE>
void Heterostructure<TYPE>::setTemperature(Temperature temp){
	for(auto layer: activeregion){
		layer.setTemperature(temp);
	}
	for(auto contact : contacts){
		contact.setTemperature(temp);
	}
	substrate->setTemperature(temp);
}


template <typename TYPE>
inline TYPE Heterostructure<TYPE>::Begin() const{
	return limitslst[0].first;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::End() const{
	return limitslst.back().second;
}


template <typename TYPE>
TYPE Heterostructure<TYPE>::getTotalWidth(){
	return limitslst.back().second-limitslst[0].first;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::gap_gamma(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->gap_gamma();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::gap_L(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->gap_L();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::gap_X(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->gap_X();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_gamma(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_gamma();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::effectiveMass_gamma() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->effectiveMass_gamma());
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_Lt(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_Lt();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::effectiveMass_Lt() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->effectiveMass_Lt());
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_Ll(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_Ll();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_Xt(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_Xt();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::effectiveMass_Xt() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->effectiveMass_Xt());
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_Xl(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_Xl();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}



template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_LDOS(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_LDOS();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_XDOS(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_XDOS();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_HHt(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_HHt();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::effectiveMass_HHt() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->effectiveMass_HHt());
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_HHl(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_HHl();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_lHt(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_lHt();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::effectiveMass_lHt() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->effectiveMass_HHt());
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::effectiveMass_lHl(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->effectiveMass_lHl();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::luttingerParam_1(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->luttingerParam_1();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::luttingerParam_2(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->luttingerParam_2();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::luttingerParam_3(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->luttingerParam_3();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::latticeparam(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->latticeparam();
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::doping(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return dopinglst[i];
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
inline TYPE Heterostructure<TYPE>::ConductionBand_gamma(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->gap_gamma()+matlst[i]->ValenceBandEnergy()+dEc[i];
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::ConductionBand_gamma() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->gap_gamma()+matlst[i]->ValenceBandEnergy()+dEc[i]);
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::ConductionBand_L(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->gap_L()+matlst[i]->ValenceBandEnergy()+dEc[i];
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::ConductionBand_L() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->gap_L()+matlst[i]->ValenceBandEnergy()+dEc[i]);
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::ConductionBand_X(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->gap_X()+matlst[i]->ValenceBandEnergy()+dEc[i];
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::ConductionBand_X() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->gap_X()+matlst[i]->ValenceBandEnergy()+dEc[i]);
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::ValenceBand_HH(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->ValenceBandEnergy()+dEv_HH[i];
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::ValenceBand_HH() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->ValenceBandEnergy()+dEv_HH[i]);
	}
	return vector;
}

template <typename TYPE>
inline TYPE Heterostructure<TYPE>::ValenceBand_lH(TYPE position) const{
	bool findpos = false;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		if(limitslst[i].first<=position){
			if(limitslst[i].second>=position){
				findpos=true;
				return matlst[i]->ValenceBandEnergy()+dEv_lH[i];
			}
		}
	}
	if(!findpos)
		throw out_of_range("The Heterostructure<TYPE> are undefined at this position");
}
template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::ValenceBand_lH() const{
	vector<TYPE> vector;
	for(unsigned int i = 0; i< limitslst.size() ;i++){
		vector.push_back(matlst[i]->ValenceBandEnergy()+dEv_lH[i]);
	}
	return vector;
}

template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::gap_gammaDiscrete(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_gamma()+dEc[currlayer]-dEv_HH[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::gap_LDiscrete(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_L()+dEc[currlayer]-dEv_HH[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::gap_XDiscrete(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_X()+dEc[currlayer]-dEv_HH[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_gamma(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_gamma();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_Lt(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_gamma();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_Ll(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_Ll();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_Xt(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_Xt();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_Xl(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_Xl();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_LDOS(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_LDOS();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_XDOS(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_XDOS();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_HHt(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_HHt();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_HHl(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_HHl();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_lHt(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_lHt();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::effectiveMass_lHl(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->effectiveMass_lHl();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::luttingerParam_1(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->luttingerParam_1();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::luttingerParam_2(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->luttingerParam_2();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::luttingerParam_3(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->luttingerParam_3();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::latticeparam(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->latticeparam();
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::doping(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=dopinglst[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::ConductionBand_gamma(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_gamma()+matlst[currlayer]->ValenceBandEnergy()
					+dEc[currlayer]+layerbeginenergy[currlayer]+Constant::e*linearelectricfield[currlayer]*(currpos-limitslst[currlayer].first)
					+electricfieldladder[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::ConductionBand_L(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_L()+matlst[currlayer]->ValenceBandEnergy()+dEc[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::ConductionBand_X(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->gap_X()+matlst[currlayer]->ValenceBandEnergy()+dEc[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}
template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::ValenceBand_HH(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->ValenceBandEnergy()+dEv_HH[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}

template <typename TYPE>
vector<TYPE> Heterostructure<TYPE>::getWidths() const{
	return widthlst;
}

template <typename TYPE>
vector<pair<TYPE,TYPE>> Heterostructure<TYPE>::getLimitsLists() const{
	return limitslst;
}

template <typename TYPE>
template <typename T>
DiscreteFunction<TYPE,T> Heterostructure<TYPE>::ValenceBand_lH(Grid1D<T>& basegrid) const{

	DiscreteFunction<TYPE,T> function = DiscreteFunction<TYPE,T>(basegrid);

	if(basegrid.getInferiorLimit()>=limitslst[0].first&&basegrid.getSuperiorLimit()<=limitslst.back().second){

		unsigned int currlayer;
		T init_pos=basegrid.getInferiorLimit();
		for(unsigned int i = 0; i< limitslst.size(); i++){
			if(limitslst[i].first<=init_pos){
				if(limitslst[i].second>=init_pos)
					currlayer = i;
			}
		}

		T currpos = init_pos;
		for(auto& itr : function){
			if(currpos>limitslst[currlayer].second)
				++currlayer;
			itr=matlst[currlayer]->ValenceBandEnergy()+dEv_lH[currlayer];
			currpos+=basegrid.getIncrement();
		}

	}
	else
		throw out_of_range("Grid outside of Heterostructure<TYPE>");

	return function;
}


static shared_ptr<Material> mixMaterial(shared_ptr<Material> one, shared_ptr<Material> two){
	if(one->getmaterialname()=="MgZnCdSe"&& two->getmaterialname()=="MgZnCdSe")
		return make_shared<MgZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()*(1-one->getCompositionQuate())+two->getComposition()*(1-two->getCompositionQuate()))/2.0,(one->getCompositionQuate()+two->getCompositionQuate())/2.0);
	else if(one->getmaterialname()=="ZnCdSe"&& two->getmaterialname()=="MgZnCdSe")
		return make_shared<MgZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()+two->getComposition()*(1-two->getCompositionQuate()))/2.0,(one->getCompositionQuate()+two->getCompositionQuate())/2.0);
	else if(one->getmaterialname()=="MgZnCdSe"&& two->getmaterialname()=="ZnCdSe")
		return make_shared<MgZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()*(1-one->getCompositionQuate())+two->getComposition())/2.0,(one->getCompositionQuate()+two->getCompositionQuate())/2.0);
	else if(one->getmaterialname()=="ZnCdSe"&& two->getmaterialname()=="ZnCdSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()+two->getComposition())/2.0);
	else if(one->getmaterialname()=="CdSe"&& two->getmaterialname()=="ZnCdSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()+two->getComposition())/2.0);
	else if(one->getmaterialname()=="ZnCdSe"&& two->getmaterialname()=="CdSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()+two->getComposition())/2.0);
	else if(one->getmaterialname()=="CdSe"&&two->getmaterialname()=="ZnSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,0.5);
	else if(one->getmaterialname()=="ZnSe"&&two->getmaterialname()=="CdSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,0.5);
	else if(one->getmaterialname()=="ZnCdSe"&& two->getmaterialname()=="ZnSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(one->getComposition()+1.0)/2.0);
	else if(one->getmaterialname()=="ZnSe"&& two->getmaterialname()=="ZnCdSe")
		return make_shared<ZnCdSe>((one->getTemperature()+two->getTemperature())/2.0,(1.0+two->getComposition())/2.0);
	else
		throw runtime_error("Mixture not valid: Roughness invalid");

}




}

#endif
