library ieee;
  use ieee.std_logic_1164.all;
  use ieee.math_real.all;

library std;       use std.textio.all;
library work;      use work.a_complex_pkg.all;

entity tb_dc_removal is
port (
    signal0DC :out complex_array_t
   
  );
end entity;

architecture sim of tb_dc_removal is
  signal longsignal  : complex_array_t;
  --signal signal0DC   : complex_array_t;
begin

  process
    file infile : text open read_mode is "C:\Users\Desktop\signal0dc.txt";
    variable L : line;
    variable re_val, im_val : real;
    variable mean_re, mean_im : real := 0.0;
    variable sum_re, sum_im   : real := 0.0;
    variable n_loaded         : integer := 0;   -- NEW: actual samples read
    variable i : integer;
  begin

    ------------------------------------------------------------------
    -- Step 1: Read up to NUM_longsignal samples (stop at EOF)
    ------------------------------------------------------------------
    for i in 0 to NUM_longsignal-1 loop
      exit when endfile(infile);                -- FIX: avoid read past EOF
      readline(infile, L);
      read(L, re_val);
      read(L, im_val);

      longsignal(i).re <= re_val;
      longsignal(i).im <= im_val;

      sum_re := sum_re + re_val;
      sum_im := sum_im + im_val;

      n_loaded := n_loaded + 1;                 -- count loaded samples
    end loop;

    wait for 0 ns;

  

    mean_re := sum_re / real(n_loaded);         -- FIX: use n_loaded
    mean_im := sum_im / real(n_loaded);

    ------------------------------------------------------------------
    -- Step 2: Subtract mean from each loaded sample
    ------------------------------------------------------------------
    for i in 0 to n_loaded-1 loop               -- FIX: loop over n_loaded
      signal0DC(i).re <= longsignal(i).re - mean_re;
      signal0DC(i).im <= longsignal(i).im - mean_im;
    end loop;

    wait;  -- end process
  end process;

end architecture;

