// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightUiTweaks from 'starlight-ui-tweaks'
import starlightScrollToTop from 'starlight-scroll-to-top'
import starlightVersions from 'starlight-versions'
import starlightSidebarTopics from 'starlight-sidebar-topics'

export default defineConfig({
  site: 'https://blinkseq.github.io/djinn',
  markdown: {
    shikiConfig: {
      themes: {
        light: 'horizon-bright',
        dark: 'horizon',
      },
    },
  },
	integrations: [
    starlight({
      customCss: ['./src/style/global.css'],
      title: '',
      logo: { src: '/src/assets/logo.svg'},
      favicon: '/favicon.png',
			social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/BLinkseq/djinn' }],
      plugins: [
        starlightUiTweaks({
        navbarLinks: [
           { label: "BLink-seq", href: "https://blinkseq.github.io"},
        ],
      }),
        starlightScrollToTop({
          showProgressRing: true,
          progressRingColor: '#b288f5',
          borderRadius: '50',
        }),
        starlightSidebarTopics([
          {
            label: 'Install',
            link: '/install',
            icon: 'download',
            items: ['install']
          },
          {
            label: 'FASTQ',
            link: '/fastq/',
            icon: 'information',
            items: [{
              label: 'FASTQ File Operations',
              items: [{ autogenerate: { directory: 'fastq' } }],
              }]
          },
          {
            label: 'SAM',
            link: '/sam/',
            icon: 'information',
            items: [{
              label: 'SAM/BAM File Operations',
              items: [{ autogenerate: { directory: 'sam' } }],
            }],
          }
        ]),
      ],
		}),
	],
});
