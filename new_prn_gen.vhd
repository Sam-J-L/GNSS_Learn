library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;
use std.textio.all;              
use IEEE.std_logic_textio.all; 
 use work.a_complex_pkg.all;

entity  new_prn_gen is
port (
    prn :out arr_2d_pr
   
  );

end entity;

architecture sim of  new_prn_gen is
	signal G1, G2 : std_logic_vector(10 downto 1) := (others => '1');
	--signal bitVal:std_logic_vector(0 to 1023*8-1);
	--signal prn : arr_2d_pr;
begin

process
  variable vG1, vG2 : std_logic_vector(10 downto 1);
  variable newG1, newG2 : std_logic;
  variable p1, p2 : std_logic;
  variable pr, i  : integer;
begin
  --report "PRN generator start";

  for pr in 0 to 31 loop                                -- pr=0..31 => PRN 1..32
    vG1 := (others => '1');                              -- seed = all ones
    vG2 := (others => '1');

    for i in 0 to 1022 loop                              -- 1023 chips
      p1 := vG1(10);                                     -- G1 output (bit 10)

      -- G2 "phase selector" taps (standard GPS PRN table 1..32)
      case pr is
        when  0 => p2 := vG2(2)  xor vG2(6);   -- PRN  1
        when  1 => p2 := vG2(3)  xor vG2(7);   -- PRN  2
        when  2 => p2 := vG2(4)  xor vG2(8);   -- PRN  3
        when  3 => p2 := vG2(5)  xor vG2(9);   -- PRN  4
        when  4 => p2 := vG2(1)  xor vG2(9);   -- PRN  5
        when  5 => p2 := vG2(2)  xor vG2(10);  -- PRN  6
        when  6 => p2 := vG2(1)  xor vG2(8);   -- PRN  7
        when  7 => p2 := vG2(2)  xor vG2(9);   -- PRN  8
        when  8 => p2 := vG2(3)  xor vG2(10);  -- PRN  9
        when  9 => p2 := vG2(2)  xor vG2(3);   -- PRN 10
        when 10 => p2 := vG2(3)  xor vG2(4);   -- PRN 11
        when 11 => p2 := vG2(5)  xor vG2(6);   -- PRN 12
        when 12 => p2 := vG2(6)  xor vG2(7);   -- PRN 13
        when 13 => p2 := vG2(7)  xor vG2(8);   -- PRN 14
        when 14 => p2 := vG2(8)  xor vG2(9);   -- PRN 15
        when 15 => p2 := vG2(9)  xor vG2(10);  -- PRN 16
        when 16 => p2 := vG2(1)  xor vG2(4);   -- PRN 17
        when 17 => p2 := vG2(2)  xor vG2(5);   -- PRN 18
        when 18 => p2 := vG2(3)  xor vG2(6);   -- PRN 19
        when 19 => p2 := vG2(4)  xor vG2(7);   -- PRN 20
        when 20 => p2 := vG2(5)  xor vG2(8);   -- PRN 21
        when 21 => p2 := vG2(6)  xor vG2(9);   -- PRN 22
        when 22 => p2 := vG2(1)  xor vG2(3);   -- PRN 23
        when 23 => p2 := vG2(4)  xor vG2(6);   -- PRN 24
        when 24 => p2 := vG2(5)  xor vG2(7);   -- PRN 25
        when 25 => p2 := vG2(6)  xor vG2(8);   -- PRN 26
        when 26 => p2 := vG2(7)  xor vG2(9);   -- PRN 27
        when 27 => p2 := vG2(8)  xor vG2(10);  -- PRN 28
        when 28 => p2 := vG2(1)  xor vG2(6);   -- PRN 29
        when 29 => p2 := vG2(2)  xor vG2(7);   -- PRN 30
        when 30 => p2 := vG2(3)  xor vG2(8);   -- PRN 31
        when 31 => p2 := vG2(4)  xor vG2(9);   -- PRN 32
        when others => p2 := '0';
      end case;

      prn(pr, i) <= p1 xor p2;                 -- C/A chip

      -- feedback: G1 taps (3,10); G2 taps (2,3,6,8,9,10)
      newG1 := vG1(3) xor vG1(10);
      newG2 := vG2(2) xor vG2(3) xor vG2(6) xor vG2(8) xor vG2(9) xor vG2(10);

      -- CORRECT shift: push right, insert new bit at bit 1 (LSB)
      vG1 := vG1(9 downto 1) & newG1;
      vG2 := vG2(9 downto 1) & newG2;
    end loop;
  end loop;
			-- report "PRN generator finish";
  wait;
end process;

end architecture;