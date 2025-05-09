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

  opts.on("-t TYPE", "--type=TYPE", "SIDIS type ('pi0' or 'pippim')") do |t|
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

unless %w[pi0 pippim].include?(options[:type])
  STDERR.puts "#{RED}Error: unsupported type '#{options[:type]}'. Supported types are 'pi0' and 'pippim'.#{RESET}"
  exit 1
end

unless File.exist?(options[:input])
  STDERR.puts "#{RED}Error: input file '#{options[:input]}' not found.#{RESET}"
  exit 1
end

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

base   = File.basename(options[:input], File.extname(options[:input]))
out    = File.join(options[:outdir], "#{base}.root")
maxev  = options[:max_events]
type   = options[:type]

# 1) Convert HIPO -> ROOT
macro = (type == 'pi0' ? 'hipo2tree_pi0' : 'hipo2tree_pippim')
hipo_cmd = %Q[clas12root -l -b -q 'macros/#{macro}.C("#{options[:input]}","#{out}",#{maxev})']
puts "#{GREEN}Running conversion:#{RESET} #{hipo_cmd}"
unless system(hipo_cmd)
  STDERR.puts "#{RED}Error: Conversion with #{macro}.C failed#{RESET}"
  exit 1
end

# 2) For pi0 only: pick GBT model and run prediction
if type == 'pi0'
  lc = options[:input].downcase
  model_type = lc.include?('outbending') ? 'outbending' : 'inbending'
  unless %w[inbending outbending].include?(model_type)
    model_type = 'inbending'
    puts "#{YELLOW}Warning: defaulting to 'inbending' model#{RESET}"
  end

  model_name = "model_rga_pass1_#{model_type}"
  model_path = File.join('src/gbt/models', model_name)
  gbt_cmd = %Q[python3 src/gbt/predict.py "#{out}" "#{model_path}" "EventTree"]
  puts "#{GREEN}Running GBT prediction with model: #{model_name}#{RESET}"
  unless system(gbt_cmd)
    STDERR.puts "#{RED}Error: GBT prediction failed#{RESET}"
    exit 1
  end
end

# 3) Run builder macro
builder = (type == 'pi0' ? 'pi0Builder' : 'pippimBuilder')
builder_cmd = %Q[clas12root -l -b -q 'macros/#{builder}.C("#{out}")']
puts "#{GREEN}Running #{builder} on:#{RESET} #{out}"
unless system(builder_cmd)
  STDERR.puts "#{RED}Error: #{builder} failed#{RESET}"
  exit 1
end

puts "#{GREEN}All steps completed successfully!#{RESET}"