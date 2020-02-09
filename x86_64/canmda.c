/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__canmda
#define _nrn_initial _nrn_initial__canmda
#define nrn_cur _nrn_cur__canmda
#define _nrn_current _nrn_current__canmda
#define nrn_jacob _nrn_jacob__canmda
#define nrn_state _nrn_state__canmda
#define _net_receive _net_receive__canmda 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define mg _p[0]
#define P _p[1]
#define i _p[2]
#define g _p[3]
#define Pca _p[4]
#define inmda _p[5]
#define gnmda _p[6]
#define iampa _p[7]
#define gampa _p[8]
#define itotal _p[9]
#define irtype _p[10]
#define f _p[11]
#define ica _p[12]
#define _g _p[13]
#define _ion_ica	*_ppvar[0]._pval
#define _ion_dicadv	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double g2_ampa;
 extern double g2_nmda;
 extern double irtype_car;
 extern double i2_ampa;
 extern double i2_nmda;
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_canmda", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
#define Area Area_canmda
 double Area = 1.11e-08;
#define k k_canmda
 double k = 1e-06;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Area_canmda", "cm2",
 "k_canmda", "mA/nA",
 "mg_canmda", "mM",
 "P_canmda", "cm/s/uS",
 "i_canmda", "mA/cm2",
 "g_canmda", "uS",
 "Pca_canmda", "cm/s",
 "inmda_canmda", "mA/cm2",
 "gnmda_canmda", "uS",
 "iampa_canmda", "mA/cm2",
 "gampa_canmda", "uS",
 "itotal_canmda", "mA/cm2",
 "irtype_canmda", "mA/cm2",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Area_canmda", &Area_canmda,
 "k_canmda", &k_canmda,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"canmda",
 "mg_canmda",
 "P_canmda",
 0,
 "i_canmda",
 "g_canmda",
 "Pca_canmda",
 "inmda_canmda",
 "gnmda_canmda",
 "iampa_canmda",
 "gampa_canmda",
 "itotal_canmda",
 "irtype_canmda",
 "f_canmda",
 0,
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	mg = 1;
 	P = 0;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _canmda_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 2);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 canmda /Users/hadasmanor/Documents/GitHubProjects/OOP_SimulationGUI/x86_64/canmda.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Ca current through NMDA receptors ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
{
 {
   P = ( 1.0 - exp ( - 65.0 * - 0.0755 ) ) / ( 10.0 * Area * 14564.0 * ( 50e-09 - ( 2e-03 * exp ( - 65.0 * - 0.0755 ) ) ) ) * k ;
   }

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = g2_nmda ;
   Pca = P * g ;
   ica = Pca * 14564.0 * v * ( 50e-09 - ( 2e-03 * exp ( v * - 0.0755 ) ) ) / ( 1.0 - exp ( v * - 0.0755 ) ) * 1.0 / ( 1.0 + ( exp ( 0.08 * - v ) * ( mg / 0.69 ) ) ) ;
   i = - Pca * 14564.0 * v * ( 50e-09 - ( 2e-03 * exp ( v * - 0.0755 ) ) ) / ( 1.0 - exp ( v * - 0.0755 ) ) * 1.0 / ( 1.0 + ( exp ( 0.08 * - v ) * ( mg / 0.69 ) ) ) ;
   gnmda = g2_nmda * 1.0 / ( 1.0 + ( exp ( 0.08 * - v ) * ( mg / 0.69 ) ) ) ;
   gampa = g2_ampa ;
   inmda = - i2_nmda ;
   iampa = - i2_ampa ;
   irtype = irtype_car ;
   itotal = i2_nmda + i2_ampa + irtype_car ;
   f = i / inmda ;
   }
 _current += ica;
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/hadasmanor/Documents/GitHubProjects/OOP_SimulationGUI/mods/canmda.mod";
static const char* nmodl_file_text = 
  "TITLE Ca current through NMDA receptors \n"
  "\n"
  ": We use this workaround mechanism to calculate the Ca current through the NMDA receptors \n"
  ": separatly from the non specific ion current through the NMDA receptors in the nmda.mod file\n"
  ": It contains:\n"
  ": \n"
  ": 1.A mechanism to caculate the Ca current through the NMDA receptor\n"
  ":   the Ca current through the NMDA receptor is added to the total Ca current \"ica(mA/cm2)\" \n"
  ":\n"
  ": 2.A balance current \"i_canmda(mA/cm2)\" (the NONSPECIFIC_CURRENT i in the \n"
  ":   code above) to the Ca current through NMDA receptors (an inward current) \n"
  ":   The balance current is needed because it has already been caculated once as a part of the \n"
  ":   total current through NMDA receptors \"i\" in the \"nmda.mod\"\n"
  ": \n"
  ": 3.Area (spine head surface area)is declared as a Global variable, and will be used in ampa.mod, nmda.mod, car.mod.\n"
  ":\n"
  ": 4.ampa, nmda and R_type current are all sent to this file as current density with the same direction of i_canmda. \n"
  ":   The itotal is just the sum of Inmda Iampa and I R_type\n"
  ":\n"
  ": Written by Lei Tian on 04/12/06 \n"
  "\n"
  "NEURON {\n"
  "	SUFFIX canmda 		\n"
  "		:will be given to the variables in this file as their family name \n"
  "	\n"
  "	USEION ca WRITE ica\n"
  "	NONSPECIFIC_CURRENT i \n"
  "	RANGE g, i, mg, inmda, gnmda, iampa, gampa, itotal, irtype, Pca, P, f\n"
  "	\n"
  "	GLOBAL Area			\n"
  "		:global varible, will be read by other files as a external one\n"
  "	\n"
  "	EXTERNAL i2_nmda, g2_nmda, i2_ampa, g2_ampa, irtype_car\n"
  "		:declare the external variables which has been declared as Global ones in nmda.mod, ampa.mod and car.mod\n"
  "	}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {                     : parameters that can be entered when function is called in cell-setup\n"
  "        dt			(ms)\n"
  "       	\n"
  "		mg   = 1	(mM)        :Mg++ concentration\n"
  "			\n"
  "		Area = 1.11e-8  (cm2)	:spine head area 1.11e-8  (cm2)\n"
  "		k = 1e-06   (mA/nA)		:transform the current from in 'nA' to in 'mA'\n"
  "\n"
  "		P           (cm/s/uS) 	:a factor to convert NMDA conductance to permeability by considering the fraction of ca current at -65mV of NMDAr is about 10% normailize it at -65mV\n"
  "}  \n"
  "\n"
  "\n"
  "ASSIGNED {	: parameters needed to solve DE\n"
  "	ica (mA/cm2)	:calcium current, which will be add to the total Ca current together with ica in 'car.mod'\n"
  "	v (mV)          :spine head membrane potential\n"
  "	i (mA/cm2)		:balance current to the ica through NMDA\n"
  "	g (uS)          :conductance of nmda(not include the effect of Mg block)\n"
  "	Pca (cm/s)		:Ca permeability of NMDA, it's obtained from gnmda by multiplied with P=0.1*gnmda*(v-e_nmda)/GHK at -65mV\n"
  "	\n"
  "	inmda (mA/cm2)	:equal to -i2_nmda which is the total nmda current's density, the direction is changed to be easier compared with i_canmda in this file.  \n"
  "	gnmda	(uS)	:cunduction of nmda(include the Mg effect), to be easier plot out by just click the 'plot what'button \n"
  "	iampa	(mA/cm2):total current of ampa, the direction is changed to be easier compared with i_canmda in this file.  \n"
  "	gampa	(uS)	:cunductance of ampa\n"
  "	itotal (mA/cm2)	:total current flow into spinehead (only the aciviated channel current is considered),the direction is chosen the same as i_canmda in this file.  \n"
  "	irtype (mA/cm2)	:r_type current\n"
  "	f               :Ca current fraction in nmda current\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "\n"
  "	P  = (1-exp(-65*-0.0755))/(10*Area*14564*(50e-09-(2e-03*exp(-65*-0.0755))))*k	:converting conductance to permaebility \n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	g = g2_nmda	:[uS]\n"
  "	Pca = P*g	:[cm/s]\n"
  "	ica = Pca*14564*v*(50e-09-(2e-03*exp(v*-0.0755)))/(1-exp(v*-0.0755))*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	:ca current density through NMDAr in [mA/cm2]\n"
  "	i = -Pca*14564*v*(50e-09-(2e-03*exp(v*-0.0755)))/(1-exp(v*-0.0755))*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	:balance current density of ca current through nmda\n"
  "	\n"
  "	:14564=(z^2*F^2)/(R*T); -0.0755 = -z*F/RT in [1/mV] where z=2,F=96500 in[C/mol], R=8.31 in[J/K*mol], T=308 in[K] \n"
  "	:and everything should be normalizied to [mV], 0.088 and 0.7474 is from our blocking experiment data fitting.\n"
  "	\n"
  "\n"
  "	gnmda=g2_nmda*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	:[uS]cunduction of nmda(include the Mg effect), to be easier plot out by just click the 'plot what'button \n"
  "	gampa=g2_ampa	:[uS]total current of ampa, the direction is changed to be easier compared with i_canmda in this file\n"
  "	inmda=-i2_nmda	:equal to -i2_nmda which is the total nmda current's density, the direction is changed to be easier compared with i_canmda in this file\n"
  "	iampa=-i2_ampa	:total current of ampa, the direction is changed to be easier compared with i_canmda in this file.  \n"
  "	\n"
  "	irtype=irtype_car	:R-type current,the direction is chosen to be easier compared with i_canmda in this file.  \n"
  "	itotal=i2_nmda+i2_ampa+irtype_car		:total current flow into spinehead (only the aciviated channel current is considered),the direction is chosen the same as i_canmda in this file. \n"
  "	f=i/inmda		:Ca current fraction in nmda current\n"
  "	\n"
  "}\n"
  "\n"
  "	\n"
  ;
#endif
