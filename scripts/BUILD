genrule(
  name="zip_outs",
  srcs=["bringup-can.sh",
        "startup.sh",
        "killall.sh",
        "//util:logger_main",
        "//can:can-dump",
        "//ui:server_main",
        "//sensor:sbus-test-run",
        "//can:can-ping",
        "//util:ping",
        "//util:pong"],
  outs=["output.zip"],
  cmd="zip -j $(location output.zip) $(SRCS)",
)
sh_binary(
  name="deploy",
  srcs=["deploy.sh"],
  data=["output.zip", "//ui:site.zip"],
  args=["$(location output.zip) $(location //ui:site.zip)"],
)
