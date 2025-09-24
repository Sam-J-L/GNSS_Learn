library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;
use IEEE.math_real.all;

library std;       use std.textio.all;
library work;      use work.a_complex_pkg.all;

entity new_a_fft_cacode is
  port(
    -- FFT outputs for all 32 PRN codes
    codeFFT : out cplx_arr_prn
  );
end entity;

architecture sim of new_a_fft_cacode is
begin
  file_read_process: process
    file infile    : text open read_mode is "C:\Users\Desktop\8192_ca.code.txt";
    variable L      : line;
    variable prnNum : integer;
    variable idx    : integer;
    variable re_v, im_v : real;
  begin
	--report"starting for loading 32 CA-code FFTs ";
    -- Read FFT data for PRNs 1..32
    for prnNum in 0 to PRN_Count-1 loop
	--if prnNum <4 then
      -- 1) skip header line "PRN X"
      readline(infile, L);

      -- 2) read MAX_SAMPLES lines of "re im"
      for idx in 0 to NFFT-1 loop
        readline(infile, L);
        read(L, re_v);
        read(L, im_v);
        -- store into the 2-D array
        codeFFT(prnNum, idx) <= (re => re_v, im => im_v);
      end loop;

      -- 3) skip the dashed separator line
      readline(infile, L);
	  --end if ;
    end loop;

    --report "done All 32 CA-code FFTs loaded." severity note;
    wait;
  end process;
end architecture;

