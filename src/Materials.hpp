/*
 * Materials.hpp
 *
 *  Created on: Jun 21, 2012
 *      Author: marcel
 */

/**
 * @file Materials.hpp
 * @author Marcel S. Claro <marcelclaro@gmail.com>
 *
 * @section License
 *
 * Copyright (c) 2012, Marcel Santos Claro
 *All rights reserved.
 *
 *Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 *1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *3. All advertising materials mentioning features or use of this software
 *   must display the following acknowledgement:
 *   This product includes software developed by the Marcel Santos Claro.
 *4. Neither the name of the Marcel Santos Claro nor the
 *   names of its contributors may be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 *THIS SOFTWARE IS PROVIDED BY Marcel Santos Claro ''AS IS'' AND ANY
 *EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL Marcel Santos Claro BE LIABLE FOR ANY
 *DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#ifndef MATERIALS_HPP_
#define MATERIALS_HPP_

#include "CMaterial.hpp"

namespace epital{


class Material;

/**
 * @class GaAs
 * @brief Galium Arsenide material
 *
 * Class for Galium Arsenide material with specific material parameters
 */
class GaAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 */
	GaAs(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/**
 * @class AlAs
 * @brief Aluminium Arsenide material
 *
 * Class for Aluminium Arsenide material with specific material parameters
 */
class AlAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 */
	AlAs(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/**
 * @class InAs
 * @brief Indium Arsenide material
 *
 * Class for Indium Arsenide material with specific material parameters
 */
class InAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 */
	InAs(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/**
 * @class AlGaAs
 * @brief Aluminium Galium Arsenide material
 *
 * Class for Aluminium Galium Arsenide material with specific material parameters
 */
class AlGaAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element (Aluminum) in alloy
	 */
	AlGaAs(Temperature temper, double composition);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/**
 * @class InGaAs
 * @brief Indium Galium Arsenide material
 *
 * Class for Indium Galium Arsenide material with specific material parameters
 */
class InGaAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy (Indium)
	 */
	InGaAs(Temperature temper, double composition);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/**
 * @class InAlAs
 * @brief Indium Aluminium Arsenide material
 *
 * Class for Indium Aluminium Arsenide material with specific material parameters
 */
class InAlAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy (Indium)
	 */
	InAlAs(Temperature temper, double composition);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};


/**
 * @class InAlGaAs
 * @brief Indium Aluminium  Gallium Arsenide material
 *
 * Class for Indium Aluminium Gallium Arsenide material with specific material parameters
 */
class InAlGaAs : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param compositionin Compositional fraction of secondary element in alloy (Indium)
	 * @param compositional Compositional fraction of secondary element in alloy (Aluminum)
	 */
	InAlGaAs(Temperature temper, double compositionin, double compositional);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/*
Missing:
class ZnCdSe;

class MgZnCdSe;

class CdSe;
*/

class ZnSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy
	 */
	ZnSe(Temperature temper);

	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class CdSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy
	 */
	CdSe(Temperature temper);

	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class ZnTe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy
	 */
	ZnTe(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class ZnSeTe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param compositionse Compositional fraction of secondary element in alloy (Selenium)
	 */
	ZnSeTe(Temperature temper, double compositionse);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class CdSeTe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param compositionte Compositional fraction of secondary element in alloy (Tellerium)
	 */
	CdSeTe(Temperature temper, double compositionte);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};



class CdTe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy
	 */
	CdTe(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class MgSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param composition Compositional fraction of secondary element in alloy
	 */
	MgSe(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class ZnCdSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param compositionzn Compositional fraction of secondary element in alloy (Zinc)
	 */
	ZnCdSe(Temperature temper, double compositionzn);

	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

class MgZnCdSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 * @param compositionzn Compositional fraction of secondary element in alloy (Zinc)
	 * @param compositionmg Compositional fraction of secondary element in alloy (Magnesium)
	 */
	MgZnCdSe(Temperature temper, double compositionzn, double compositionmg);


	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

inline double ValenceBandOffset(const Material& mat_a,const Material& mat_b){
	return mat_b.ValenceBandEnergy()-mat_a.ValenceBandEnergy();
}


/**
 * @class GaSe
 * @brief Galium Selenide material
 *
 * Class for Galium Selenide material with specific material parameters
 */
class GaSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 */
	GaSe(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};

/**
 * @class InSe
 * @brief Indium Selenide material
 *
 * Class for Indium Selenide material with specific material parameters
 */
class InSe : public Material{
public:
	/**
	 * Initialize all materials parameters
	 * @param temper Initial temperature
	 */
	InSe(Temperature temper);
	/**
	 * Change temperature and refresh material parameters with new values
	 * @param temper New temperature
	 */
	void setTemperature(Temperature temper);
};





}

#endif /* MATERIALS_HPP_ */
