// Subscribe after startup so instant navigation also typesets new page content.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  },
  startup: {
    ready() {
      MathJax.startup.defaultReady();
      MathJax.startup.promise.then(() => {
        if (typeof document$ === "undefined") return;
        let pending = Promise.resolve();
        document$.subscribe(() => {
          pending = pending.then(() => {
            MathJax.typesetClear();
            MathJax.texReset();
            return MathJax.typesetPromise();
          }).catch(error => console.error("MathJax typesetting failed:", error));
        });
      });
    }
  }
};
