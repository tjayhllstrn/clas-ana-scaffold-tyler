#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'

# ANSI color codes
RESET  = "\e[0m"
GREEN  = "\e[1;32m"
YELLOW = "\e[1;33m"
RED    = "\e[1;31m"

options = { max_events: 100 }

parser = OptionParser.new do |opts|
  opts.banner = "#{GREEN}Usage: run_single_hipo.rb -t TYPE -i FILE -o DIR [options]#{RESET}"

  opts.on("-t TYPE", "--type=TYPE", "SIDIS type (only 'pi0' supported)") do |t|
    options[:type] = t
  end

  opts.on("-i FILE", "--input=FILE", "Input .hipo file") do |f|
    options[:input] = f
  end

  opts.on("-o DIR", "--outdir=DIR", "Output directory") do |d|
    options[:outdir] = d
  end

  opts.on("-n N", "--max-events=N", Integer,
          "Max events (default #{options[:max_events]}, use -1 for all)") do |n|
    options[:max_events] = n
  end

  opts.on("-h", "--help", "Show this help") do
    puts opts
    exit
  end
end

begin
  parser.parse!
  missing = %i[type input outdir].select { |k| options[k].nil? }
  unless missing.empty?
    STDERR.puts "#{RED}Missing options: #{missing.join(', ')}#{RESET}"
    STDERR.puts parser
    exit 1
  end
rescue OptionParser::ParseError => e
  STDERR.puts "#{RED}#{e.message}#{RESET}"
  STDERR.puts parser
  exit 1
end

unless options[:type] == 'pi0'
  STDERR.puts "#{RED}Error: unsupported type '#{options[:type]}'. Only 'pi0' is supported.#{RESET}"
  exit 1
end

unless File.exist?(options[:input])
  STDERR.puts "#{RED}Error: input file '#{options[:input]}' not found.#{RESET}"
  exit 1
end

# --- ensure output directory exists or offer to create it ---
unless Dir.exist?(options[:outdir])
  print "#{YELLOW}Output directory '#{options[:outdir]}' does not exist. Create it? [y/N]: #{RESET}"
  ans = STDIN.gets.chomp.downcase
  if %w[y yes].include?(ans)
    FileUtils.mkdir_p(options[:outdir])
    puts "#{GREEN}Created directory #{options[:outdir]}#{RESET}"
  else
    STDERR.puts "#{RED}Aborted: output directory not found.#{RESET}"
    exit 1
  end
end

base  = File.basename(options[:input], File.extname(options[:input]))
out   = File.join(options[:outdir], "#{base}.root")
maxev = options[:max_events]

# 1) Convert HIPO → ROOT
hipo_cmd = %Q[clas12root -l -b -q "macros/hipo2tree_pi0.C(\\"#{options[:input]}\\",\\"#{out}\\",#{maxev})"]
puts "#{GREEN}Running conversion:#{RESET} #{hipo_cmd}"
unless system(hipo_cmd)
  STDERR.puts "#{RED}Error: Failed to run hipo2tree_pi0.C#{RESET}"
  exit 1
end

# 2) Pick GBT model based on "inbending"/"outbending" in path
lc = options[:input].downcase
if lc.include?('outbending')
  model_type = 'outbending'
elsif lc.include?('inbending')
  model_type = 'inbending'
else
  model_type = 'inbending'
  puts "#{YELLOW}Warning: input path does not contain 'inbending' or 'outbending'; defaulting to 'inbending'.#{RESET}"
end

model_name = "model_rga_pass1_#{model_type}"
model_path = File.join('src/gbt/models', model_name)

# 3) Run GBT prediction
gbt_cmd = %Q[python3 src/gbt/predict.py "#{out}" "#{model_path}" "EventTree"]
puts "#{GREEN}Running GBT prediction with model: #{model_name}#{RESET}"
unless system(gbt_cmd)
  STDERR.puts "#{RED}Error: GBT prediction failed#{RESET}"
  exit 1
end

# 4) Run pi0Builder on the same ROOT file
pi0_cmd = %Q[clas12root -l -b -q "macros/pi0Builder.C(\\"#{out}\\")"]
puts "#{GREEN}Running pi0Builder on:#{RESET} #{out}"
unless system(pi0_cmd)
  STDERR.puts "#{RED}Error: pi0Builder failed#{RESET}"
  exit 1
end

puts "#{GREEN}All steps completed successfully!#{RESET}"
