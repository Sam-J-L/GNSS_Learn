-- File: a_test.vhd
library ieee;
  use ieee.std_logic_1164.all;
  use ieee.math_real.all;

library std;       use std.textio.all;
library work;      use work.a_complex_pkg.all;

entity fft_try is
 -- port (
  --  IQ1_out : out cplx_arr_2d;  -- [FB][k], k=0..8183
  --  IQ2_out : out cplx_arr_2d
  --);
end entity;

architecture sim of fft_try is

  signal IQ1_out :cplx_arr_2d;
  signal IQ2_out :cplx_arr_2d;
  signal codeFFT_sig :cplx_arr_prn;
  signal conv1_sig     :cplx_arr_3d;
  signal conv2_sig     :cplx_arr_3d;
  signal corr1_time : real_arr_3d;  -- [pr][fb][n]
  signal corr2_time,results : real_arr_3d;
  signal max_results: real_arr_2d;
  signal peakSize,fftMax: real_arr_1d;
  signal frequencyBinIndex,fftMaxIndex:integer_arr_1d;
  signal codePhaseRange :integer_arr ;
  signal resultPeakSize :other_arr ;
  signal secondPeakSize : real_arr_prn;
 
  signal prn1,caCode: arr_2d_pr;
  signal signal10DC: complex_array_t;
  signal codeValueIndex : int_array_t(0 to arrayLength - 1);
   signal tsignal10DC,xCarrier : complex_array_tt ;
  signal longCaCode     : bit_array_t(0 to arrayLength - 1);
  --signal fftxc  : real_vec(0 to FFT_NUMPTS-1);
--=============just try===================================
	-- Constants (static => OK for declarations)
	constant ARRAY_LEN  : natural := xCarrier'length;  -- e.g., 10*NFFT
	constant FFT_LOG2   : natural := ceil_log2(ARRAY_LEN) + 3;  -- +3 == ×8
	--constant FFT_NUMPTS : natural := 2**FFT_LOG2; 
	constant fftNumPts  :integer  :=8*2**(ceil_log2(xCarrier'length));
 constant FS    : real := 8.192e6;	-- 8 * 2^ceil(log2(L))
 constant REV_KEEP : natural := 2;
 constant uniqFftPts:integer := ceil_vhdl( (real(fftNumPts) + 1.0) / 2.0 );
    signal fftxc  : real_vec(0 to fftNumPts-1);
	signal fftFreqBins : real_vec1(0 to uniqFftPts-1);
	signal fftFreqBinsRev,tem : real_vec1(0 to (uniqFftPts-3));
	-- Zero-padded buffer, FFT output, and |FFT|
 signal acqResults_peakMetric,acqResults_carrFreq: real_arr_prn;
 signal acqResults_codePhase:integer_arr_1d;
--==================================

  --------------------------------------------------------------------------

 
begin
 ----------------------
 fft_ca : entity work.new_a_fft_cacode
    
    port map (
      codeFFT => codeFFT_sig
    );
prngen: entity work.new_prn_gen
	port map (
		prn => prn1
		);
sognalDC: entity work.tb_dc_removal
	port map (
		signal0DC => signal10DC
		);
---------------------------
  main_proc: process
    -- files
    file f1       : text open read_mode is "C:\Users\KAR\Desktop\new_sample_signal1.txt";	-- will change to differant form 
    file f2       : text open read_mode is "C:\Users\KAR\Desktop\new_sample_signal2.txt";
    variable L        : line;
    variable re_tmp   : real;
    variable im_tmp   : real;

    -- buffers & grids
    variable phaseP   : real_arr;
    variable frqCand  : frq_arr;
    variable sig1v    : cplx_arr_sig;
    variable sig2v    : cplx_arr_sig;
    variable I1_arr   : real_arr;
    variable Q1_arr   : real_arr;
    variable I2_arr   : real_arr;
    variable Q2_arr   : real_arr;

    -- FFT outputs for each segment (first 8184 bins)
    variable FFT1_bins : cplx_arr_sig;
    variable FFT2_bins : cplx_arr_sig;

    -- indices/temps
    variable cnt1, cnt2  : integer := 0;
    variable fb, idx     : integer;
    variable angle, cre, cim : real;
	-- convoulation 
	variable a_re, a_im   : real;
    variable c_re, c_im   : real;
    variable r_re, r_im   : real;
	---------------------
	variable bins_line : cplx_arr_sig;
    variable td_line,out_corr1,out_corr2,res   : real_arr;
	variable m_idx,frq_idx,t : integer;
	variable mag : real;
	variable acqresult1,acqresult2 : real;
	variable max_res : real;
	variable rangeIndex,indx2,src: integer;
	--===================
	variable I_arr, Q_arr : real_arr;        -- 0..NFFT-1
	variable Xbins        : cplx_arr_sig;    -- FFT output
	--variable fftxc        : real_arr;        -- |FFT|, same meaning as MATLAB fftxc
	--==============================================
	variable x_pad  : cplx_vec(0 to fftNumPts-1);
	variable X_fft  : cplx_vec(0 to fftNumPts-1);
	--==============================================
	variable fftMax: real; 
	variable fftMaxIndex,x : integer;
	 variable w : integer ;
	 
	variable samplesPerCodeChip,excludeRangeIndex1,excludeRangeIndex2: integer;
  begin
    --report"process is starting";
	--report ""&integer'image(fftNumPts);
    ----------------------------------------------------------------
    -- 1) phasePoints = n*2π·TS
    ----------------------------------------------------------------
    for idx in 0 to NFFT-1 loop
      phaseP(idx) := real(idx) * TWOPI * TS;
    end loop;

    ----------------------------------------------------------------
    -- 2) read 1-ms segments
    ----------------------------------------------------------------
    while (cnt1 < NFFT) and not endfile(f1) loop
      readline(f1, L); read(L, re_tmp); read(L, im_tmp);
      sig1v(cnt1).re := re_tmp;  sig1v(cnt1).im := im_tmp;
      cnt1 := cnt1 + 1;
    end loop;

    while (cnt2 < NFFT) and not endfile(f2) loop
      readline(f2, L); read(L, re_tmp); read(L, im_tmp);
      sig2v(cnt2).re := re_tmp;  sig2v(cnt2).im := im_tmp;
      cnt2 := cnt2 + 1;
    end loop;

    ----------------------------------------------------------------
    -- 3) per-Doppler: mix to baseband, then FFT (replaces DFT loops)
    ----------------------------------------------------------------
	
 for pr in 0 to PRN_Count-1 loop    
	for fb in 0 to FRQ_BINS-1 loop
		  frqCand(fb) := FREQ_START + real(fb) * FREQ_STEP;

--========================= generate I/Q signals & FFT to I/Q & CONV ==========================
		  for idx in 0 to NFFT-1 loop
			angle       := frqCand(fb) * phaseP(idx);
			cre         := cos(angle);
			cim         := sin(angle);

			I1_arr(idx) := cre*sig1v(idx).re - cim*sig1v(idx).im;
			Q1_arr(idx) := cre*sig1v(idx).im + cim*sig1v(idx).re;

			I2_arr(idx) := cre*sig2v(idx).re - cim*sig2v(idx).im;
			Q2_arr(idx) := cre*sig2v(idx).im + cim*sig2v(idx).re;
		
		  end loop;
		  --for idx in 0 to 5 loop
		--		report"I1_arr("&integer'image(idx)&")= "&real'image(I1_arr(idx));
		--end loop;
		  -- (b) FFT on segment 1   -- FFT: replaces DFT on seg1
		  fft8192_from_IQ(I1_arr, Q1_arr, FFT1_bins);
		  for idx in 0 to NFFT-1 loop
			IQ1_out(fb, idx) <= (re => FFT1_bins(idx).re, im => FFT1_bins(idx).im);
		  end loop;
			wait for 0 ns ;
		  -- (c) FFT on segment 2   -- FFT: replaces DFT on seg2
		  fft8192_from_IQ(I2_arr, Q2_arr, FFT2_bins);
		  for idx in 0 to NFFT-1 loop
			IQ2_out(fb, idx) <= (re => FFT2_bins(idx).re, im => FFT2_bins(idx).im);
		  end loop;
		  wait for 0 ns; 
	for idx in 0 to NFFT-1 loop
        -- unpack IQ1_sig for this PRN and bin
        a_re := IQ1_out(fb, idx).re;  
        a_im := IQ1_out(fb, idx).im;   

        -- unpack the PRN FFT for this PRN
        c_re := codeFFT_sig(pr, idx).re;  
        c_im := codeFFT_sig(pr, idx).im;  

        -- complex multiply → conv1_sig(pr,fb,idx)
        r_re := (a_re * c_re) - (a_im * c_im);
        r_im := (a_re * c_im) + (a_im * c_re);
        conv1_sig(pr, fb, idx) <= (re => r_re, im => r_im);

        -- same for IQ2_sig → conv2_sig(pr,fb,idx)
        a_re := IQ2_out(fb, idx).re;   
        a_im := IQ2_out(fb, idx).im;   
        r_re := (a_re * c_re) - (a_im * c_im);
        r_im := (a_re * c_im) + (a_im * c_re);
        conv2_sig(pr, fb, idx) <= (re => r_re, im => r_im);  
      end loop;  -- idx
	  --===========================================================
			  wait for 0 ns;

--=========================FFT of conv1_sig(pr, fb, :)=====================================
		for k in 0 to NFFT-1 loop
		  bins_line(k) := conv1_sig(pr, fb, k);
		end loop;
		ifft8192(bins_line, td_line);
		for k in 0 to NFFT-1 loop
		  corr1_time(pr, fb, k) <= td_line(k);
		end loop;

		-- IFFT of conv2_sig(pr, fb, :)
		for k in 0 to NFFT-1 loop
		  bins_line(k) := conv2_sig(pr, fb, k);
		end loop;
		ifft8192(bins_line, td_line);
		for k in 0 to NFFT-1 loop
		  corr2_time(pr, fb, k) <= td_line(k);
		end loop;
		wait for 0 ns; 
--====================================CONV Done ==============================

		for k in 0 to NFFT-1 loop
		  out_corr1(k) := corr1_time(pr, fb, k);
		  out_corr2(k) := corr2_time(pr, fb, k);
		end loop;
		-----------------------
		max(out_corr1,mag,m_idx);
		acqresult1:= mag;
		max(out_corr2,mag,m_idx);
		acqresult2:= mag;
		-----------------------
		for n in 0 to NFFT-1 loop
			if acqresult1 > acqresult2 then 
				results(pr,fb,n) <= out_corr1(n);
			else
				results(pr,fb,n) <= out_corr2(n);
			end if;	
	    end loop;
		--============================================
		wait for 0 ns ;
		for n in 0 to NFFT-1 loop 
		res(n):=results(pr,fb,n);
			max(res,mag,m_idx);
		end loop;
		--wait for 0 ns;
		max_results(pr,fb) <= mag;
		wait for 0 ns ;
--===================== frequency bins & peak size ===========================
		for n in 0 to FRQ_BINS-1 loop
		res(n):= max_results(pr,n);
		end loop; 
		max(res,mag,frq_idx);
		--report "freq=  "& integer'image(frq_idx);
		peakSize(pr)<=mag;
	   frequencyBinIndex(pr)<=frq_idx;
--===================== code phase =========================================================
	for n in 0 to NFFT-1 loop 
		res(n):=results(pr,frq_idx,n);
	--	max(res,mag,m_idx);
	end loop;
	max(res,mag,m_idx);
	acqResults_codePhase(pr)<=m_idx;
	wait for 0 ns ; 
--===================== exclude range around the peak ===========================
	samplesPerCodeChip := (NFFT / codeFreqBasis1);
	excludeRangeIndex1 := acqResults_codePhase(pr) - samplesPerCodeChip;
    excludeRangeIndex2 := acqResults_codePhase(pr) + samplesPerCodeChip;
	----------------------------------------------------------------------------
	rangeIndex:=0;
	if excludeRangeIndex1 < 2 then 
	for i in excludeRangeIndex2 to NFFT loop 
		codePhaseRange(rangeIndex)<= i;
		rangeIndex:= rangeIndex+1;
	end loop ;
	
	for i in 0 to NFFT + excludeRangeIndex1 loop 
		codePhaseRange(rangeIndex)<= i;
		rangeIndex:= rangeIndex+1;
	end loop;
	--------------
	elsif excludeRangeIndex2 >= NFFT then 
	rangeIndex:=0;
	for i in excludeRangeIndex2-NFFT to excludeRangeIndex1 loop 
	codePhaseRange(rangeIndex)<= i;
		rangeIndex:= rangeIndex+1;
	end loop ;
	------------------------------------------------------------------
	else
	rangeIndex:=0;
	for i in 0 to excludeRangeIndex1 loop
		codePhaseRange(rangeIndex)<= i;
		rangeIndex:= rangeIndex+1;
	end loop;
	for i in excludeRangeIndex2 to NFFT-1 loop
		codePhaseRange(rangeIndex)<= i;
		rangeIndex:= rangeIndex+1;
	end loop;
	end if ;
	wait for 0 ns ;
--===================Find the second highest correlation peak in the same freq. bin==========================================
	for i in 0 to 8176 loop
	resultPeakSize(i) <= results(pr,frq_idx,codePhaseRange(i));
	end loop;
	wait for  0 ns;
	max(resultPeakSize,mag,m_idx);
	secondPeakSize(pr)<=mag;
	wait for 0 ns;
	--indx2:=m_idx;
--=================== % If the result is above threshold, then there is a signal ...===========================
 acqResults_peakMetric(pr)<= peakSize(pr)/secondPeakSize(pr);
	
	if  acqResults_peakMetric(pr) > acqThreshold then 
				--report "there is a signal at PRN = "&integer'image(pr);
			for i in 0 to 1022 loop 
				caCode(pr,i)<=prn1(pr,i);
			end loop;
			
			for i in 1 to arrayLength  loop
				codeValueIndex(i-1) <= floor_vhdl((ts * real(i)) / (1.0 / codeFreqBasis));
			end loop;
				wait for 0 ns ; 
			
		-----------longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
			for i in 0 to arrayLength - 1 loop
				longCaCode(i) <= caCode(pr, codeValueIndex(i) mod 1023);
			end loop;
				wait for 0 ns ;
--               Remove C/A code modulation from the original signal -(Using detected C/A code phase)-----xCarrier = ...signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ... .* longCaCode
			
			for k in 0 to arrayLength-1 loop
			  -- wrap-around into source buffer
			  src := signal10DC'low + ((acqResults_codePhase(pr) - signal10DC'low + k) mod signal10DC'length); --0+((6234-0+0)mod 44....)=6

			  -- copy one complex sample into the k-th position of the window
			  tsignal10DC(k) <= signal10DC(src);
			end loop;
				wait for 0 ns ;
			for k in 0 to arrayLength-1 loop
				if longCaCode(k) = '1' then
			 
					xCarrier(k) <= tsignal10DC(k);
				else
			 
				xCarrier(k) <= (re => -tsignal10DC(k).re,
							  im => -tsignal10DC(k).im);
				
				end if;
			end loop;
			wait for 0 ns ;
			--wait for 0 ns ;
			------------Find the next highest power of two and increase by 8x 
			--fftNumPts := 8*(2^(nextpow2(length(xCarrier))));
			--fftNumPts:=8*2**(ceil_log2(xCarrier'length));
				--report"fftNumPts= "&integer'image(fftNumPts);


		-- Zero-pad xCarrier to length FFT_NUMPTS
		for k in 0 to fftNumPts-1 loop
			 if k <= xCarrier'high then
			x_pad(k) := xCarrier(k);
			 else
			x_pad(k).re := 0.0;
			x_pad(k).im := 0.0;
			 end if;
		end loop;

			 --=====--FFT: X_fft = fft(x_pad)
		fft_pow2(x_pad, X_fft);
--
			--Magnitude: fftxc = abs(X_fft)
		for k in 0 to fftNumPts-1 loop
			fftxc(k) <= sqrt( X_fft(k).re*X_fft(k).re + X_fft(k).im*X_fft(k).im );
		end loop;
		 max(fftxc,mag,frq_idx);
		 fftMax:=mag;
	     fftMaxIndex:=frq_idx;
		 wait for 0 ns ;
		--uniqFftPts := ceil_vhdl( (real(fftNumPts) + 1.0) / 2.0 ); define above 
		--=============================================
		for k in 0 to uniqFftPts-1 loop
		
		fftFreqBins(k)<= real(K)*(FS/real(fftNumPts));
		
		end loop;
			wait for 0 ns ;
			if fftMaxIndex > uniqFftPts then 
				w:=0;
				if (fftNumPts rem 2) = 0 then             -- even length: exclude DC and Nyquist
				  for k in (uniqFftPts-2) downto 1 loop
					exit when w > fftFreqBinsRev'high;
					fftFreqBinsRev(w) <= -fftFreqBins(k);
					w := w + 1;
					
				  end loop;
				  -- [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
				  w:=0;
				  for k in uniqFftPts to fftxc'length-1 loop
					exit when w > tem'high;
					tem(w)<=fftxc(k);
					w:= w+1;
				  end loop;
				   max(tem,mag,frq_idx);
					fftMax:=mag;
					fftMaxIndex:=frq_idx;
					acqResults_carrFreq(pr)<=-fftFreqBinsRev(fftMaxIndex);
				else                                     -- odd length: skip DC
				  for k in (uniqFftPts-2) downto 1 loop
					exit when w > fftFreqBinsRev'high;
					fftFreqBinsRev(w) <= -fftFreqBins(k);
					w := w + 1;                
				  end loop;
				   w:=1;
				  for k in uniqFftPts+3 to fftxc'length loop
					exit when w > tem'high;
					tem(w)<=fftxc(k);
					w:= w+1;
				  end loop;
				   max(tem,mag,frq_idx);
					fftMax:=mag;
					fftMaxIndex:=frq_idx;
					acqResults_carrFreq(pr)<=-fftFreqBinsRev(fftMaxIndex);
					
				end if ;
			end if ;
	--report "acqResults_carrFreq"&real'image(acqResults_carrFreq(pr));
		--report"x= "&integer'image(uniqFftPts);
		 --report"fftMax= "&real'image(mag);
		 --report"fftMaxIndex= "&integer'image(frq_idx);
		wait for 0 ns ;		
		--=================================================
    else
	acqResults_carrFreq(pr)<=0.0;
	end if;						--		 if (peakSize/secondPeakSize) > settings.acqThreshold


	
	end loop; --FB
	t:=pr;
	report"pr = "&integer'image(t)&" done";
   end loop; --pr
   -- report "IQ1_out/IQ2_out complete (FFT-based)." severity note;
    wait;
  end process;
end architecture;
--if pr=0 and fb=0 then 
		--report"acqresult1= "&real'image(acqresult1)&"  acqresult2= "&real'image(acqresult2);
		--end if ;
		--wait;
	--====================
--report "samplesPerCodeChip=  "& integer'image(samplesPerCodeChip);
	--report "excludeRangeIndex1=  "& integer'image(excludeRangeIndex1);
	--report "excludeRangeIndex2=  "& integer'image(excludeRangeIndex2);
