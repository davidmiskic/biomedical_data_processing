function [idx] = QRSDetect(fileName, m, normCnst)
%function QRSDetect(fileName, m, normCnst)
  S = load(fileName);
  sig0 = S.val(1,:);
  sig1 = S.val(2,:);
  sigLen = size(sig1,2);
  CONST_SAMPLES_PER_S = 250;
  
  % CONDITIONER
  % first filter, b is for x, a is for y. put minus for values for y on
  % right and minus for values for x on left
  b1 = [1 2 1];
  a1 = 4;
  % second filter 
  b2 = [1 -(2*cos((60*pi)/125)) 1];
  a2 = 1;
  % third filter
  b3 = [1 0 0 0 0 0 -1];
  a3 = 1;
  % energy collector filter
  bc = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
  ac = 20;

  sig0 = filter(b1, a1, sig0);
  sig0 = filter(b2, a2, sig0);
  sig0 = filter(b3, a3, sig0);
  sig0 = abs(sig0);
  sig1 = filter(b1, a1, sig1);
  sig1 = filter(b2, a2, sig1);
  sig1 = filter(b3, a3, sig1);
  sig1 = abs(sig1);
  
  sig = (sig0+sig1)/2;
  sig_secondary = filter(bc, ac, sig.^2);

  % SECONDARY DETECTOR
  ET = mean(sig_secondary(1:CONST_SAMPLES_PER_S));  % init threshold
  i = 1;
  idx_secondary = [];
  last_r = 0; % where last peak was located
  rr_interval = (CONST_SAMPLES_PER_S/1000)*200; % at least this much samples between QRSs
  amplitude_means = []; % store last 8 peaks
  width_limit = [(CONST_SAMPLES_PER_S/1000)*16 (CONST_SAMPLES_PER_S/1000)*500]; % width of event between 16 and 500ms
  after_detection_lower = (CONST_SAMPLES_PER_S/1000)*200; % 200ms after event lower threshold
  
  while i <= sigLen
      if sig_secondary(i) > ET
          end_event = i+1;
          peak_val = sig_secondary(i);
          peak_i = i;
          while end_event <= sigLen && sig_secondary(end_event) >= ET  % event end only after drop below ET
              if sig_secondary(end_event) > peak_val
                  peak_val = sig_secondary(end_event);
                  peak_i = end_event;
              end
              end_event = end_event + 1;
          end
          event_width = end_event - i;
          if (last_r == 0 || (peak_i - last_r > rr_interval)) && ... % more than 200ms from previous QRS
              (event_width >= width_limit(1) && event_width <= width_limit(2)) && ... % event between 16 and 500 ms long
              (size(amplitude_means, 2) == 0 || (mean(amplitude_means)*0.1 <= peak_val && mean(amplitude_means)*6 >= peak_val)) % between 10% and 600% of mean of 8 previous
              %ET = 1*(ET + 0.5*peak_val);
              ET = 2*(ET + 0.25*peak_val);
              %ET = 4*(ET + 0.5*peak_val);
              last_r = peak_i;
              idx_secondary = [idx_secondary peak_i];
              if size(amplitude_means, 2) > 8
                  amplitude_means = [amplitude_means(1:7) peak_val];
              else 
                  amplitude_means = [amplitude_means peak_val];
              end
          end
          i = end_event;
      elseif size(idx_secondary, 2) ~= 0 && i - idx_secondary(end) == after_detection_lower % 200ms after detection decrease threshold
          ET = 0.2 * ET;
      elseif size(idx_secondary, 2) ~= 0 && i - idx_secondary(end) == CONST_SAMPLES_PER_S % if 1 second with no detection decrease threshold
          ET = 0.5 * ET;
          if size(amplitude_means, 2) ~= 0 && amplitude_means(end) > mean(amplitude_means) % if last mean lower, decrease last 8 mean, else increase
            amplitude_means = amplitude_means.*1.25;
          elseif size(amplitude_means, 2) ~= 0 && amplitude_means(end) < mean(amplitude_means)
              amplitude_means = amplitude_means.*0.75;
          end
      else 
          i = i + 1;
      end
  end
  disp(size(idx_secondary, 2) + " SEC detect");  

  % MAIN DETECTOR
  BT = 100;
  DT = 100;
  time_interval = 180; % in ms
  samples_interval = round((CONST_SAMPLES_PER_S / 1000) * time_interval);
  i = 1;
  idx = [];
  while i <= sigLen && samples_interval > 0
      if sig(i) > DT % possible event
          n_cross = 0;
          peak_val = sig(i);
          peak_i = i;
          j = i+1;
          while j <= sigLen && j <= i+samples_interval % count crossings in interval
            if ((sig(j) > DT && sig(j-1) <= DT) || (sig(j-1) > DT && sig(j) <= DT)) && j-1 ~= i
            %if (sig(j) > DT && sig(j-1) <= DT) && j-1 ~= i
                n_cross = n_cross + 1;
            end
            if sig(j) > peak_val
                peak_val = sig(j);
                peak_i = j;
            end
            j = j + 1;
          end
          %disp(n_cross);
          if 2 <= n_cross && n_cross <= 4 % QRS detected
              %disp(idx);
              idx = [idx peak_i];
              BT = 0.75 * BT + 0.25 * peak_val;
              DT = max([BT 0.5*peak_val]);
              %i = i+1;
              i = j;
          elseif n_cross == 0 % nothing
              i = i + 1;
              BT = 0.5*BT;
              DT = BT;
          elseif n_cross > 4 % noise
              i = i + 1;
              BT = BT * 1.5;
              DT = max([BT 0.5*peak_val]);
          elseif n_cross == 1 % undecided
              % did secondary find anything between i and j
              contains = false;
              for possible_i = i:j
                  if any(idx_secondary == possible_i)
                      contains = true;
                      break;
                  end
              end
              if contains
                  idx = [idx peak_i];
                  %i = i+1;
                  i = j;
              else
                  i = i + 1;
              end
              DT = 0.75 * DT;
              DT = max([DT 0.5*BT]);
          else
              i = i + 1;
          end
          BT = min([1200 BT]);
          DT = min([1200 DT]);
      else
          i = i + 1;
      end
  end % main detector 
  disp(size(idx, 2) + " PRIM detect"); 
  
  for i=1:size(idx_secondary, 2)
      elem = idx_secondary(i);
      contains = false;
      for possible_i = elem-50:elem+50
          if any(idx == possible_i)
              contains = true;
              break;
          end
      end
      if contains == false
          idx = [idx elem];
      end
  end
  
  disp(size(idx, 2) + " FINAL"); 
end
