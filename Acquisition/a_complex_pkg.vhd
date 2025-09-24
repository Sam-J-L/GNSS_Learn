-- File: a_complex_pkg.vhd
library ieee;
library work; 
library std;   
  use ieee.std_logic_1164.all;
  use ieee.numeric_std.all;
  use ieee.math_real.all;

package a_complex_pkg is
  constant MAX_SAMPLES   : integer := 8184;   -- one C/A code period
  constant FRQ_BINS      : integer := 29;     -- Doppler bins
  constant NFFT          : integer := 8192;   -- DFT length
  constant codeFreqBasis : real := 1.023e6;
  constant codeFreqBasis1 : integer := 1023;
  constant acqThreshold  : real := 2.5;
  constant PRN_Count     : integer:=7;
  constant arrayLength    : integer := 10 * NFFT;
  -------------------------
  constant IF_FREQ     : real    := 38400.0;       -- Hz
  constant SEARCH_BAND : real    := 14.0;          -- kHz
  constant FREQ_START  : real    := IF_FREQ - (SEARCH_BAND/2.0)*1000.0;
  constant FREQ_STEP   : real    := 500.0;         -- 0.5 kHz steps
  ----------------
    -- Sampling/time constants
  constant FS    : real := 8.192e6;       -- as you defined
  constant TS    : real := 1.0/FS;
  constant TWOPI : real := 2.0 * 3.141592653589793;
  constant NUM_longsignal : integer := 442368;

  -- one complex sample
  type cplx_t is record
    re : real;
    im : real;
  end record;
  --
  

  -- I/Q buffer (1 ms segment): 8184 samples
  type cplx_arr_sig is array(0 to NFFT-1) of cplx_t;
  type complex_array_t is array (0 to NUM_longsignal-1) of cplx_t;
   type complex_array_tt is array (0 to arrayLength-1) of cplx_t;
  
  ---------------------------------------------------

  -- NFFT-length array for DFTs & code FFT
  type cplx_arr     is array(0 to NFFT-1)       of cplx_t;
  type cplx_idx		is array(0 to 3) of cplx_t;

  -- phasePoints or I/Q‐magnitude arrays
  --type real_arr     is array(0 to NFFT-1) of real;
  type integer_arr     is array(0 to 8176) of integer; -- i think should determine the size auto after cal but in mine time do that 
 -- type other_arr     is array(0 to 8176) of real;

  -- Doppler frequency grid
  type frq_arr      is array(0 to FRQ_BINS-1)   of real;
  ---------------------
   type cplx_arr_2d     is array(0 to FRQ_BINS-1, 0 to NFFT-1) of cplx_t;
 --------------------------------------
 type cplx_arr_22d     is array(0 to FRQ_BINS-1, 0 to NFFT-1) of real;
 type arr_2d_pr        is array(0 to 31,0 to 1022) of std_logic;
 -------------------------------------------------
  type max_arr     is array(0 to 28) of real;
  type max_arr_multi_prns     is array(0 to PRN_Count-1, 0 to 28) of real;
  ---------------------------------------------------------
  --for new codes 
  type cplx_arr_3d    is array(0 to PRN_Count-1,       0 to FRQ_BINS-1,  0 to NFFT-1) of cplx_t;
  type real_arr_3d    is array(0 to PRN_Count-1,       0 to FRQ_BINS-1,  0 to NFFT-1) of real;
  type real_arr_2d    is array(0 to PRN_Count-1,       0 to FRQ_BINS-1) of real;
  type real_arr_1d    is array(0 to PRN_Count-1) of real;
  type integer_arr_1d    is array(0 to PRN_Count-1) of integer;
  type real_arr_prn    is array(0 to PRN_Count-1) of real;
  type cplx_arr_prn is array(0 to PRN_Count-1, 0 to NFFT-1) of cplx_t;
  --===========================subtype=================================
  type real_vec is array (integer range <>) of real;

-- Subtypes (examples)
subtype real_arr        is real_vec(0 to NFFT-1);   -- 8192
subtype other_arr       is real_vec(0 to 8176);     -- 8169
subtype frq_line_arr    is real_vec(0 to FRQ_BINS-1);
subtype real_vec1 		is real_vec; 
  --------------------any length--------------------------------
  type int_array_t is array(natural range <>) of integer;
  type bit_array_t is array(natural range <>) of std_logic;
  type cplx_vec is array (natural range <>) of cplx_t;
  --type real_vec1 is array (natural range <>) of real;
   G_PRNS : natural := PRN_Count;
    G_NFFT : natural := NFFT;  
  type rom2d_t is array (0 to G_PRNS-1, 0 to G_NFFT-1) of coeff_t;

  --==========================FFT & IFFT ==============================
  function round_real(x : real) return integer;
  function floor_vhdl(x : real) return integer;
  function ceil_vhdl(x : real) return integer;
  function ceil_log2(n : positive) return natural;
  function bitrev(i : integer; bits : integer) return integer;
   procedure fft8192_from_IQ(
    I_arr : in  real_arr;
    Q_arr : in  real_arr;
    bins  : out cplx_arr_sig
  );
  procedure ifft8192(
  Xbins : in  cplx_arr_sig;   -- length NFFT spectrum
  mag  : out real_arr   -- length NFFT time-domain result
);
procedure max(
  Xbins : in  real_vec;
  max_value  : out real;
  m_idx : out integer
) ;
procedure fft_pow2(
  X_in  : in  cplx_vec;   
  Y_out : out cplx_vec    
  );
end package;

--==========================================
package body a_complex_pkg is
 -- Implementation of the round_real function
  function round_real(x : real) return integer is
  begin
    if x >= 0.0 then
      return integer(x + 0.5);
    else
      return integer(x - 0.5);
    end if;
  end function round_real;
  --======floor function 
  function floor_vhdl(x : real) return integer is
  variable i : integer;
begin
  i := integer(x);
  if real(i) > x then
    return i - 1;
  else
    return i;
  end if;
end function;
--==========ceil function=============================
function ceil_vhdl(x : real) return integer is
  variable i : integer := integer(x);  -- truncates toward 0
begin
  if x > real(i) then
    return i + 1;
  else
    return i;
  end if;
end function;
--==smallest k s.t. 2**k >= n-----------
function ceil_log2(n : positive) return natural is
	  variable v : natural := n - 1;
	  variable k : natural := 0;
	begin
	  while v > 0 loop
		v := v / 2;
		k := k + 1;
	  end loop;
	  return k;
end function;
  --==================================================
  -- bit-reverse an integer of 'bits' width (for NFFT=8192, bits=13)
  function bitrev(i : integer; bits : integer) return integer is
    variable v : integer := i;
    variable r : integer := 0;
  begin
    for b in 1 to bits loop
      r := (r * 2) + (v mod 2);
      v := v / 2;
    end loop;
    return r;
  end function;
  --====================================================
  
  --====================================================
    -- FFT procedure: pack I/Q → zero-pad → in-place radix-2 DIT FFT → copy first 8184 bins
  procedure fft8192_from_IQ(
    I_arr  : in  real_arr;           -- 0..8183
    Q_arr  : in  real_arr;           -- 0..8183
    bins   : out cplx_arr_sig        -- 0..8183 (first MAX_SAMPLES bins of FFT)
  ) is
    constant LOG2N : integer := 13;  -- log2(8192)
    variable X     : cplx_arr;       -- 8192 work buffer
    variable tmp_re, tmp_im : real;
    variable m, half_m, kblk, j     : integer;
    variable angle, wr, wi, tr, ti  : real;
    variable ur, ui                 : real;
	variable r : integer;
  begin
    -- 1) pack and zero-pad
    for n in 0 to NFFT-1 loop
      X(n).re := I_arr(n);
      X(n).im := Q_arr(n);
    end loop;
   -- for n in MAX_SAMPLES to NFFT-1 loop
   --   X(n).re := 0.0;
   --   X(n).im := 0.0;
   -- end loop;

    -- 2) bit-reversal permutation
    for n in 0 to NFFT-1 loop
       r:= bitrev(n, LOG2N);
      if r > n then
        tmp_re := X(n).re;
		tmp_im := X(n).im;
        X(n).re := X(r).re; 
		X(n).im := X(r).im;
        X(r).re := tmp_re; 
        X(r).im := tmp_im;
      end if;
    end loop;

    -- 3) iterative radix-2 DIT FFT (DFT sign: e^{-j 2π kn/N})
    for s in 1 to LOG2N loop
      m      := 2**s; --2^s=2---M=4
      half_m := m/2; -- =1,,=2
      -- blocks of size m
      kblk := 0; --kblk = block start index in the radix-2 DIT FFT.
      while kblk < NFFT loop
        for j in 0 to half_m-1 loop
          angle := -TWOPI * real(j) / real(m);  -- twiddle W_m^j
          wr    := cos(angle);--1 
          wi    := sin(angle);--0

          -- t = W * X[k+j+half]
          tr := wr*X(kblk+j+half_m).re - wi*X(kblk+j+half_m).im;
          ti := wr*X(kblk+j+half_m).im + wi*X(kblk+j+half_m).re;

          -- u = X[k+j]
          ur := X(kblk+j).re;--x(0)
          ui := X(kblk+j).im;--x(0)

          -- butterflies
          X(kblk+j).re          := ur + tr; -- x(0)=x(0)+x(1) -- x(2)
          X(kblk+j).im          := ui + ti; -- x(0) -- x(2)=x(2)+x(3)
          X(kblk+j+half_m).re   := ur - tr; -- x(1) -- x(3)
          X(kblk+j+half_m).im   := ui - ti; -- x(1) -- x(3)
        end loop;
        kblk := kblk + m;
      end loop;
    end loop;
	--report "kblk= " & integer'image(kblk);
    -- 4) copy first MAX_SAMPLES bins to output
    for k in 0 to NFFT-1 loop
      bins(k).re := X(k).re;
      bins(k).im := X(k).im;
    end loop;
	--report"done fft";
  end procedure;
  --==============================================================
  procedure ifft8192(
  Xbins : in  cplx_arr_sig;
  mag  : out real_arr
) is
  variable tmpI : real_arr;         -- re{conj(X)} = re{X}
  variable tmpQ : real_arr;         -- im{conj(X)} = -im{X}
  variable Y    : cplx_arr_sig;     -- FFT{conj(X)}
  variable xout :  cplx_arr_sig;
  constant INVN : real := 1.0 / real(NFFT);
begin
  -- conj the input spectrum
  for k in 0 to NFFT-1 loop
    tmpI(k) := Xbins(k).re;
    tmpQ(k) := -Xbins(k).im;
  end loop;

  -- forward FFT on conj(X)
  fft8192_from_IQ(tmpI, tmpQ, Y);

  -- x[n] = (1/N) * conj( Y[n] )
  for n in 0 to NFFT-1 loop
    xout(n).re :=  Y(n).re * INVN;
    xout(n).im := -Y(n).im * INVN;
	mag(n):=xout(n).re*xout(n).re+xout(n).im*xout(n).im;
  end loop;
  --report"done ifft";
end procedure;
--========================max===================================
procedure max(
  Xbins     : in  real_vec;
  max_value : out real;
  m_idx     : out integer
) is
  variable mv : real;
  variable mi : integer;
begin
  assert Xbins'length > 0 ;

  mv := Xbins(Xbins'left);
  mi := Xbins'left;

  for n in Xbins'range loop
    if Xbins(n) > mv then
      mv := Xbins(n);
      mi := n;
    end if;
  end loop;

  max_value := mv;
  m_idx     := mi;
end procedure;
----------===================new fft===========================================================
procedure fft_pow2(
  X_in  : in  cplx_vec;   -- any length N=2^L, any natural range (e.g., 0..N-1)
  Y_out : out cplx_vec    -- same range/length as X_in
) is
  constant N      : natural := X_in'length;
  constant L2     : natural := ceil_log2(N);
  constant BASE   : natural := X_in'low;

  -- Work buffer with same bounds as X_in
  variable X : cplx_vec(X_in'range);

  -- temps
  variable tmp_re, tmp_im : real;
  variable ur, ui, tr, ti : real;
  variable angle, wr, wi   : real;
  variable m, half_m       : natural;
  variable kblk, j         : natural;
  variable n0, n1          : natural;  -- absolute indices
  variable r               : natural;  -- bit-reversed index (0-based)
begin
  --Sanity checks
-- assert (2**L2 = N)
    --report "fft_pow2: length is not a power of two"
  -- severity failure;
 --assert (Y_out'length = N)
   -- report "fft_pow2: X_in and Y_out lengths differ"
   --severity failure;

  -- 1) copy input to work buffer
  for i in X'range loop
    X(i) := X_in(i);
  end loop;

  -- 2) bit-reversal on 0-based positions
for n in 0 to N-1 loop
  r := bitrev(n, integer(L2));
  if r > n then
    n0 := n;         -- <— use BASE here
    n1 := r;         -- <— and here
    tmp_re := X(n0).re;  tmp_im := X(n0).im;
    X(n0).re := X(n1).re;  X(n0).im := X(n1).im;
    X(n1).re := tmp_re;    X(n1).im := tmp_im;
  end if;
end loop;

  -- 3) iterative radix-2 DIT FFT (DFT sign: e^{-j 2π kn/N})
  for s in 1 to L2 loop
    m      := 2**s;
    half_m := m/2;
    kblk   := 0;
    while kblk < N loop
      for j in 0 to half_m-1 loop
        angle := -TWOPI * real(j) / real(m);
        wr    := cos(angle);
        wi    := sin(angle);

        n0 := BASE + kblk + j;
        n1 := n0 + half_m;

        -- t = W * X[n1]
        tr := wr*X(n1).re - wi*X(n1).im;
        ti := wr*X(n1).im + wi*X(n1).re;

        -- u = X[n0]
        ur := X(n0).re;
        ui := X(n0).im;

        -- butterflies
        X(n0).re := ur + tr;
        X(n0).im := ui + ti;
        X(n1).re := ur - tr;
        X(n1).im := ui - ti;
      end loop;
      kblk := kblk + m;
    end loop;
  end loop;

  -- 4) write to output
  for i in Y_out'range loop
    Y_out(i) := X(i);
  end loop;
end procedure;
--==============================================================
------------------------------------------------------
end package body;
