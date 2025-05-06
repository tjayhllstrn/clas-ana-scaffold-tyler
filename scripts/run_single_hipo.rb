#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'

options = {
  max_events: 100
}

parser = OptionParser.new do |opts|
  opts.banner = "Usage: run_single_hipo.rb -t TYPE -i FILE -o DIR [options]"

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
  missing = %i[type input outdir].select { |key| options[key].nil? }
  unless missing.empty?
    STDERR.puts "Missing options: #{missing.join(', ')}"
    STDERR.puts parser
    exit 1
  end
rescue OptionParser::ParseError => e
  STDERR.puts e
  STDERR.puts parser
  exit 1
end

unless options[:type] == 'pi0'
  STDERR.puts "Error: unsupported type '#{options[:type]}'. Only 'pi0' is supported."
  exit 1
end

unless File.exist?(options[:input])
  STDERR.puts "Error: input file '#{options[:input]}' not found."
  exit 1
end

# --- ensure output directory exists or offer to create it ---
unless Dir.exist?(options[:outdir])
  print "Output directory '#{options[:outdir]}' does not exist. Create it? [y/N]: "
  ans = STDIN.gets.chomp.downcase
  if ans == 'y' || ans == 'yes'
    FileUtils.mkdir_p(options[:outdir])
    puts "Created directory #{options[:outdir]}"
  else
    abort("Aborted: output directory not found.")
  end
end

# --- decide overwrite vs append if the output file already exists ---
base  = File.basename(options[:input], File.extname(options[:input]))
out   = File.join(options[:outdir], "#{base}.root")
maxev = options[:max_events]


cmd = %Q[clas12root -l -b -q "macros/hipo2tree_pi0.C(\\"#{options[:input]}\\",\\"#{out}\\",#{maxev})"]
puts "Running: #{cmd}"
system(cmd) or abort("Failed to run hipo2tree_pi0.C")
